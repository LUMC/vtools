"""
vtools.gcoverage
~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
import itertools
from typing import Generator, Iterable, List, NamedTuple, Tuple, Union

import cyvcf2

import numpy as np


class Region(NamedTuple):
    chr: str
    start: int
    end: int

    def __str__(self):
        return "{0}:{1}-{2}".format(self.chr, self.start, self.end)


def coverage_for_gvcf_record(record: cyvcf2.Variant, maxlen: int = 15000
                             ) -> List[int]:
    """
    Get coverage for gvcf record per base

    Some records may be huge, especially those around centromeres.
    Therefore, there is a maxlen argument. Maximally `maxlen` values
    are returned.
    """
    size = record.end - record.start
    try:
        dp = record.format("DP")[0][0]
    except TypeError:
        dp = 0
    if size <= maxlen:
        return [dp] * size
    else:
        return [dp] * maxlen


def gq_for_gvcf_record(record: cyvcf2.Variant, maxlen: int = 15000
                       ) -> List[int]:
    """
    Some records may be huge, especially those around centromeres.
    Therefore, there is a maxlen argument. Maximally `maxlen` values
    are returned.
    """
    size = record.end - record.start
    gq = record.gt_quals[0]
    if size <= maxlen:
        return [gq] * size
    else:
        return [gq]*maxlen


def qualmean(quals: np.ndarray) -> float:
    """
    Credit:
    https://gigabaseorgigabyte.wordpress.com/2017/06/26/averaging-basecall-quality-scores-the-right-way/
    https://git.lumc.nl/klinische-genetica/capture-lumc/vtools/issues/3
    """
    return -10*np.log10(np.mean(np.power(10, quals/-10)))


def fractions_at_least(values: np.ndarray,
                       boundaries: Iterable[Union[int, float]]
                       ) -> List[float]:
    total = values.size
    bins = [0] + list(boundaries) + [np.iinfo(values.dtype).max]
    counts, _ = np.histogram(values, bins)
    # Example 10,20,30,40. Will give us bins 0-9,10-19,20-29,30-39,40-max.
    # 40-max contains the count of everything at least 40. Everything at least
    # 30 should include bins 30-39 and 40-max. Etc. The results we seek are
    # therefore cumulative from high to low. Hence the below function.
    reverse_cumulative_counts = list(itertools.accumulate(reversed(counts)))
    cumulative_fractions = [count / total for count in
                            reversed(reverse_cumulative_counts)]
    # We are not interested in the lowest bin. It will be 1.000 anyway as it
    # contains all samples.
    return cumulative_fractions[1:]


class CovStats(NamedTuple):
    median_dp: float
    median_gq: float
    mean_dp: float
    mean_gq: float
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
        perc_at_least_dp = [fraction * 100 for fraction in
                            fractions_at_least(coverages,
                                               (10, 20, 30, 50, 100))]
        perc_at_least_gq = [fraction * 100 for fraction in
                            fractions_at_least(gq_qualities,
                                               (10, 20, 30, 50, 90))]
        return cls(np.median(coverages),
                   np.median(gq_qualities),
                   np.mean(coverages),
                   np.mean(gq_qualities),
                   *perc_at_least_dp,
                   *perc_at_least_gq)

    @classmethod
    def header(cls):
        return "\t".join(cls.__annotations__.keys())

    def __str__(self):
        # 8 width, 2 decimals for nice right aligned values
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
        if len(contents) < 11:
            raise ValueError("refFlat line must have at least 11 fields")
        gene = contents[0]
        transcript = contents[1]
        contig = contents[2]
        if "-" in contents[3].strip():
            forward = False
        else:
            forward = True
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
        regs = []
        for s, e in zip(self.exon_starts, self.exon_ends):
            regs.append(Region(self.contig, s, e))
        return regs

    @property
    def cds_exons(self) -> List[Tuple[int, Region]]:
        regs = []
        for i, (s, e) in enumerate(zip(self.exon_starts, self.exon_ends)):
            if s < self.cds_start:
                s = self.cds_start
            if e > self.cds_end:
                e = self.cds_end
            reg = Region(self.contig, s, e)
            if reg.end <= reg.start:   # utr exons
                continue
            regs.append((i, reg))
        return regs


def file_to_refflat_records(filename: str) -> Generator[RefRecord, None, None]:
    with open(filename, "rt") as file_h:
        for line in file_h:
            yield RefRecord.from_line(line)


def feature_to_vcf_records(feature: List[Region], sample_vcfs: List[cyvcf2.VCF]
                           ) -> Generator[cyvcf2.Variant, None, None]:
    for sample_vcf in sample_vcfs:
        for region in feature:
            for record in sample_vcf(str(region)):
                yield record


def gvcf_records_to_coverage_array(gvcf_records: Iterable[cyvcf2.Variant],
                                   maxlen: int = 15000
                                   ) -> np.ndarray:
    coverages = itertools.chain(coverage_for_gvcf_record(gvcf_record, maxlen)
                                for gvcf_record in gvcf_records)
    # Coverages higher than 65535 will be incorrectly recorded
    # but this will save memory and these occurrences should be
    # very rare. Since coverage is never below 0 use an unsigned
    # integer.
    return np.fromiter(coverages, dtype=np.uint16)


def gvcf_records_to_gq_array(gvcf_records: Iterable[cyvcf2.Variant],
                             maxlen: int = 15000
                             ) -> np.ndarray:
    coverages = itertools.chain(gq_for_gvcf_record(gvcf_record, maxlen)
                                for gvcf_record in gvcf_records)
    # GQ can never be higher than 99 and not lower than 0.
    # uint8 with values 0-255 is appropriate. It is an integer in
    # GVCF file format.
    return np.fromiter(coverages, dtype=np.uint8)


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
                coverage = gvcf_records_to_coverage_array(records)
                gq_quals = gvcf_records_to_gq_array(records)
                covstats = CovStats.from_coverages_and_gq_qualities(
                    coverage, gq_quals)
                exon = i + 1 if refflat_record.forward else total_exons - i
                yield f"{gene}\t{transcript}\t{exon}\t{str(covstats)}"
    else:
        yield "gene\ttranscript\t" + CovStats.header()
        for refflat_record in file_to_refflat_records(refflat_file):
            regions = [x[1] for x in refflat_record.cds_exons]
            records = feature_to_vcf_records(regions, gvcf_readers)
            coverage = gvcf_records_to_coverage_array(records)
            gq_quals = gvcf_records_to_gq_array(records)
            covstats = CovStats.from_coverages_and_gq_qualities(
                    coverage, gq_quals)
            yield (f"{refflat_record.gene}\t{refflat_record.transcript}\t"
                   f"{str(covstats)}")
