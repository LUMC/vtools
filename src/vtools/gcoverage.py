# MIT License
#
# Copyright (c) 2018, 2020 Leiden University Medical Center
# Copyright (c) 2018 Sander Bollen
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Calculate coverage statistics using GVCF files.
"""

import argparse
import os
from typing import Generator, Iterable, List, NamedTuple, Tuple, Union

import cyvcf2  # type: ignore

import numpy as np


class Region(NamedTuple):
    """Represents a VCF region"""
    chr: str
    start: int
    end: int

    def __str__(self):
        return "{0}:{1}-{2}".format(self.chr, self.start, self.end)

    def __len__(self):
        return self.end - self.start + 1  # +1 since the position is one-based.


def qualmean(quals: np.ndarray) -> float:
    """
    Credit:
    https://gigabaseorgigabyte.wordpress.com/2017/06/26/averaging-basecall-quality-scores-the-right-way/
    https://git.lumc.nl/klinische-genetica/capture-lumc/vtools/issues/3
    """
    return -10 * np.log10(np.mean(np.power(10, quals / -10)))


class CovStats(NamedTuple):
    """Class representing a line in a CovStats TSV."""
    mean_dp: float
    mean_gq: float
    median_dp: float
    median_gq: float
    perc_at_least_10_dp: float
    perc_at_least_20_dp: float
    perc_at_least_30_dp: float
    perc_at_least_50_dp: float
    perc_at_least_100_dp: float
    perc_at_least_10_gq: float
    perc_at_least_20_gq: float
    perc_at_least_30_gq: float
    perc_at_least_50_gq: float
    perc_at_least_90_gq: float

    @classmethod
    def from_coverages_and_gq_qualities(cls,
                                        coverage_list: List[int],
                                        gq_quality_list: List[float]):
        """Generate a CovStats object from an array of coverages and genome
        qualities."""
        total_coverages = len(coverage_list)
        total_qualities = len(gq_quality_list)

        # np.fromiter is faster than np.array in this case. Also specifying the
        # length allows allocating the array at once in memory, which is
        # faster. Use int64 and float64 which are compatible with python
        # integers and floats. Using smaller 32-bit types does not improve
        # performance noticably.
        coverages = np.fromiter(
            coverage_list, dtype=np.int64, count=total_coverages)
        gq_qualities = np.fromiter(
            gq_quality_list, dtype=np.float64, count=total_qualities)

        # np.count_nonzero checks the __bool__ method of objects. (In python2
        # this was called __nonzero__). So this is the appropriate way for
        # counting booleans. (Sum also works but is less correct and slower).
        perc_at_least_dp = [
            np.count_nonzero(np.greater_equal(coverages, at_least))
            / total_coverages * 100  # percentage calculation.
            for at_least in (10, 20, 30, 50, 100)]
        perc_at_least_gq = [
            np.count_nonzero(np.greater_equal(gq_qualities, at_least))
            / total_qualities * 100
            for at_least in (10, 20, 30, 50, 90)]

        # Explicit float conversion needed as numpy does not guarantee a float.
        return cls(float(np.mean(coverages)),
                   qualmean(gq_qualities),
                   float(np.median(coverages)),
                   float(np.median(gq_qualities)),
                   *perc_at_least_dp,
                   *perc_at_least_gq)

    @classmethod
    def header(cls, compact: bool = False):
        if compact:
            return ("mean_dp\tmean_gq\tmedian_dp\tmedian_gq\t%dp>=10\t"
                    "%dp>=20\t%dp>=30\t%dp>=50\t%dp>=100\t%gq>=10\t%gq>=20\t"
                    "%gq>=30\t%gq>=50\t%gq>=90")
        return "\t".join(cls.__annotations__.keys())

    def __str__(self):
        # 2 decimals is enough precision.
        return "\t".join("{:.2f}".format(value) for value in self)


class RefRecord(NamedTuple):
    gene: str
    transcript: str
    contig: str
    start: int
    end: int
    cds_start: int
    cds_end: int
    exon_starts: List[int]
    exon_ends: List[int]
    forward: bool

    @classmethod
    def from_line(cls, line):
        contents = line.strip().split("\t")
        if len(contents) != 11:
            raise ValueError(f"refFlat line must have exactly 11 fields. "
                             f"Error on line: {line}.")
        gene = contents[0]
        transcript = contents[1]
        contig = contents[2]
        strand = contents[3].strip()
        if strand == "+":
            forward = True
        elif strand == "-":
            forward = False
        else:
            raise ValueError(f"Invalid strand: '{strand}'. Strand should be "
                             f"'+' or '-'. Error on line: {line}. ")
        start = int(contents[4])
        end = int(contents[5])
        cds_start = int(contents[6])
        cds_end = int(contents[7])
        exon_starts = [int(x) for x in contents[9].split(",")[:-1]]
        exon_ends = [int(x) for x in contents[10].split(",")[:-1]]
        return cls(gene, transcript, contig, start, end, cds_start,
                   cds_end, exon_starts, exon_ends, forward)

    @property
    def exons(self) -> List[Region]:
        return [Region(self.contig, s, e)
                for s, e in zip(self.exon_starts, self.exon_ends)]

    @property
    def cds_exons(self) -> List[Region]:
        regs = []
        for s, e in zip(self.exon_starts, self.exon_ends):
            if s < self.cds_start:
                s = self.cds_start
            if e > self.cds_end:
                e = self.cds_end
            reg = Region(self.contig, s, e)
            # Since positioning is one-based reg.end == reg.start is a valid
            # region (of one basepair). Hence < instead of <=.
            # Coding sequences that are smaller than 3 basepairs are probably
            # non-existent in the wild, but it happens in the test data.
            if reg.end < reg.start:  # utr exons
                continue
            regs.append(reg)
        return regs


def file_to_refflat_records(filename: Union[str, os.PathLike]
                            ) -> Generator[RefRecord, None, None]:
    with open(filename, "rt") as file_h:
        for line in file_h:
            yield RefRecord.from_line(line)


def feature_to_coverage_and_quality_lists(feature: List[Region],
                                          vcfs: List[cyvcf2.VCF]
                                          ) -> Tuple[List[int], List[float]]:
    """
    Wrapper around region_and_vcf_to_coverage_and_quality_lists that returns
    the values for multiple regions in multiple vcfs.
    """
    depths: List[int] = []
    gen_quals: List[float] = []
    for vcf in vcfs:
        for region in feature:
            covs, quals = region_and_vcf_to_coverage_and_quality_lists(
                region, vcf)
            depths.extend(covs)
            gen_quals.extend(quals)
    return depths, gen_quals


def region_and_vcf_to_coverage_and_quality_lists(
        region: Region,
        vcf: cyvcf2.VCF
) -> Tuple[List[int], List[float]]:
    """
    Gets the coverages and qualities for a particular region in a VCF.
    """
    depths: List[int] = []
    gen_quals: List[float] = []
    for variant in vcf(str(region)):
        # max and min to make sure the positions outside of the region are
        # not considered.
        # Technically max and min should only be considered for the first
        # and last record respectively. But doing it for every record does not
        # affect the outcome while also being correct for the edge case with
        # only one record (which is both first and last) and using
        # significantly less code. Also the extra logic needed is probably
        # slower than running max and min on everything.
        start = max(region.start - 1, variant.start)  # -1 for one-based pos.
        end = min(region.end, variant.end)
        size = end - start
        gq = variant.gt_quals[0]
        try:
            dp = variant.format("DP")[0][0]
        except TypeError:
            # If there is no DP field, None will be returned.
            # None[0] returns TypeError: 'NoneType' object is not subscriptable
            dp = 0
        depths.extend([dp] * size)
        gen_quals.extend([gq] * size)
    return depths, gen_quals


def refflat_record_to_regions(refflat_record: RefRecord,
                              region_of_interest: str) -> List[Region]:
    if region_of_interest == "transcript":
        return [Region(refflat_record.contig, refflat_record.start,
                       refflat_record.end)]
    elif region_of_interest == "transcript_cds_exons":
        return refflat_record.cds_exons
    else:
        raise ValueError(f"Unsupported region: {region_of_interest}")


def refflat_and_gvcfs_to_tsv(refflat_file: str,
                             gvcfs: Iterable[str],
                             region_of_interest: str,
                             compact_header: bool = False,
                             ) -> Generator[str, None, None]:
    gvcf_readers = [cyvcf2.VCF(gvcf) for gvcf in gvcfs]
    per_exon = region_of_interest == "exon"
    if per_exon:
        yield "gene\ttranscript\texon\t" + CovStats.header(compact_header)
        for refflat_record in file_to_refflat_records(refflat_file):
            total_exons = len(refflat_record.exons)
            gene = refflat_record.gene
            transcript = refflat_record.transcript
            for i, region in enumerate(refflat_record.exons):
                coverage, gq_quals = (
                    feature_to_coverage_and_quality_lists(
                        [region], gvcf_readers))
                covstats = CovStats.from_coverages_and_gq_qualities(
                    coverage, gq_quals)
                exon = i + 1 if refflat_record.forward else total_exons - i
                yield f"{gene}\t{transcript}\t{exon}\t{str(covstats)}"
    else:
        yield "gene\ttranscript\t" + CovStats.header(compact_header)
        for refflat_record in file_to_refflat_records(refflat_file):
            regions = refflat_record_to_regions(refflat_record,
                                                region_of_interest)
            coverage, gq_quals = (
                feature_to_coverage_and_quality_lists(
                    regions, gvcf_readers))
            if coverage and gq_quals:
                # If cds exons are processed, sometimes coverage and gq_quals
                # will be empty lists as some transcripts do not have protein-
                # coding exons. In that case skip.
                covstats = CovStats.from_coverages_and_gq_qualities(
                    coverage, gq_quals)
                yield (f"{refflat_record.gene}\t{refflat_record.transcript}\t"
                       f"{str(covstats)}")


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate coverage statistics using one or more "
                    "single-sample gvcf files. The statistics are calculated "
                    "over the intervals in the refflat file(s).")
    parser.add_argument("input_gvcf", metavar="GVCF", nargs="+", type=str,
                        help="The GVCF file(s) over which the coverage is "
                             "calculated. When multiple files are given the "
                             "average is calculated.  Multisample gvcf files "
                             "are not supported.")
    refflat_grp = parser.add_mutually_exclusive_group(required=True)
    refflat_grp.add_argument("-R", "--refflat-file", type=str,
                             help="Path to refFlat file.")
    refflat_grp.add_argument("-Z", "--refflat-zip", type=str,
                             help="Zip file containing multiple refflats. The "
                                  "output will be a zipfile with multiple TSV "
                                  "files.")
    method_grp = parser.add_mutually_exclusive_group()
    parser.set_defaults(region_of_interest='exon')
    method_grp.add_argument("--per-exon", action="store_const",
                            dest="region_of_interest", const='exon',
                            help="Collect metrics per exon")
    method_grp.add_argument("--per-transcript", action="store_const",
                            dest="region_of_interest", const='transcript',
                            help="Collect metrics per transcript")
    method_grp.add_argument("--per-transcript-cds-exons", action="store_const",
                            dest="region_of_interest",
                            const='transcript_cds_exons',
                            help="Collect metrics per transcript, only "
                                 "considering the exons in the coding region.")
    parser.add_argument("-s", "--short-column-names", action="store_true",
                        help="Print shorter column names for easier viewing "
                             "on a terminal.")
    return parser


def main():
    args = argument_parser().parse_args()
    for line in refflat_and_gvcfs_to_tsv(args.refflat_file, args.input_gvcf,
                                         args.region_of_interest,
                                         args.short_column_names):
        print(line)
