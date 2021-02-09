"""
vtools.gcoverage
~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
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


def qualmean(quals: np.ndarray) -> float:
    """
    Credit:
    https://gigabaseorgigabyte.wordpress.com/2017/06/26/averaging-basecall-quality-scores-the-right-way/
    https://git.lumc.nl/klinische-genetica/capture-lumc/vtools/issues/3
    """
    return -10*np.log10(np.mean(np.power(10, quals/-10)))


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
                                        coverages: np.ndarray,
                                        gq_qualities: np.ndarray):
        """Generate a CovStats object from an array of coverages and genome
        qualities."""
        total_coverages = coverages.size
        total_qualities = gq_qualities.size
        # numpy.greater_equal returns a list of Booleans. But since true==1
        # these can be summed for the total count.
        perc_at_least_dp = [
            np.greater_equal(coverages, at_least).sum() / total_coverages
            * 100 for at_least in (10, 20, 30, 50, 100)]
        perc_at_least_gq = [
            np.greater_equal(gq_qualities, at_least).sum() / total_qualities
            * 100 for at_least in (10, 20, 30, 50, 90)]

        # Explicit float conversion needed as numpy does not guarantee a float.
        return cls(float(np.mean(coverages)),
                   qualmean(gq_qualities),
                   float(np.median(coverages)),
                   float(np.median(gq_qualities)),
                   *perc_at_least_dp,
                   *perc_at_least_gq)

    @classmethod
    def header(cls):
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
            if reg.end <= reg.start:   # utr exons
                continue
            regs.append(reg)
        return regs


def file_to_refflat_records(filename: Union[str, os.PathLike]
                            ) -> Generator[RefRecord, None, None]:
    with open(filename, "rt") as file_h:
        for line in file_h:
            yield RefRecord.from_line(line)


def feature_to_vcf_records(feature: List[Region], sample_vcfs: List[cyvcf2.VCF]
                           ) -> Generator[cyvcf2.Variant, None, None]:
    for sample_vcf in sample_vcfs:
        for region in feature:
            for record in sample_vcf(str(region)):
                yield record


def gvcf_records_to_coverage_and_quality_arrays(
        gvcf_records: Iterable[cyvcf2.Variant],
        maxlen: int = 15000
        ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get coverage and genome quality for all gvcf records per base

    Some records may be huge, especially those around centromeres.
    Therefore, there is a maxlen argument. Maximally `maxlen` values
    are used per variant.
    """
    depths: List[int] = []
    gen_quals: List[int] = []
    for record in gvcf_records:
        size = record.end - record.start
        try:
            dp = record.format("DP")[0][0]
        except TypeError:
            dp = 0
        gq = record.gt_quals[0]
        record_size = min(size, maxlen)  # Limit size to maxlen
        depths.extend([dp] * record_size)
        gen_quals.extend([gq] * record_size)
    # np.fromiter is faster than np.array in this case. Also specifying the
    # length allows allocating the array at once in memory, which is faster.
    # Use int64 which is compatible with python integer. Using smaller
    # integers does not improve performance noticably.
    return (np.fromiter(depths, dtype=np.int64, count=len(depths)),
            np.fromiter(gen_quals, dtype=np.int64, count=len(gen_quals)))


def refflat_and_gvcfs_to_tsv(refflat_file: str,
                             gvcfs: Iterable[str],
                             per_exon=False
                             ) -> Generator[str, None, None]:
    gvcf_readers = [cyvcf2.VCF(gvcf) for gvcf in gvcfs]
    if per_exon:
        yield "gene\ttranscript\texon\t" + CovStats.header()
        for refflat_record in file_to_refflat_records(refflat_file):
            total_exons = len(refflat_record.exons)
            gene = refflat_record.gene
            transcript = refflat_record.transcript
            for i, region in enumerate(refflat_record.exons):
                records = feature_to_vcf_records([region], gvcf_readers)
                coverage, gq_quals = (
                    gvcf_records_to_coverage_and_quality_arrays(records))
                covstats = CovStats.from_coverages_and_gq_qualities(
                    coverage, gq_quals)
                exon = i + 1 if refflat_record.forward else total_exons - i
                yield f"{gene}\t{transcript}\t{exon}\t{str(covstats)}"
    else:
        yield "gene\ttranscript\t" + CovStats.header()
        for refflat_record in file_to_refflat_records(refflat_file):
            records = feature_to_vcf_records(
                refflat_record.cds_exons, gvcf_readers)
            coverage, gq_quals = (
                gvcf_records_to_coverage_and_quality_arrays(records))
            covstats = CovStats.from_coverages_and_gq_qualities(
                    coverage, gq_quals)
            yield (f"{refflat_record.gene}\t{refflat_record.transcript}\t"
                   f"{str(covstats)}")
