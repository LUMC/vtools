"""
vtools.gcoverage
~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
import itertools
from itertools import chain
from typing import List, Optional, Tuple, NamedTuple, Iterable, Union, \
    Generator

import cyvcf2

import numpy as np


class Region(NamedTuple):
    chr: str
    start: int
    end: int

    def __str__(self):
        return "{0}:{1}-{2}".format(self.chr, self.start, self.end)


def coverage_for_gvcf_record(record: cyvcf2.Variant, maxlen: int = 15000) -> List[int]:
    """
    Get coverage for gvcf record per base

    Some records may be huge, especially those around centromeres.
    Therefore, there is a maxlen argument. Maximally `maxlen` values
    are returned.
    """
    start = record.start
    end = record.end
    size = end - start
    try:
        dp = record.format("DP")[0][0]
    except TypeError:
        dp = 0
    if size < maxlen:
        return [dp]*(end-start)
    else:
        return [dp]*maxlen


def gq_for_gvcf_record(record: cyvcf2.Variant, maxlen: int = 15000) -> List[int]:
    """
    Some records may be huge, especially those around centromeres.
    Therefore, there is a maxlen argument. Maximally `maxlen` values
    are returned.
    """
    start = record.start
    end = record.end
    size = end - start
    gq = record.gt_quals[0]
    if size < maxlen:
        return [gq]*(end-start)
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
    # We are not interested in the lowest bin. It will be 1.000_ anyway.
    return cumulative_fractions[1:]


class CovStats(object):
    def __init__(self, records):
        self.records = records
        self.__coverages = None
        self.__gq_qualities = None

    @property
    def coverages(self) -> np.ndarray:
        if self.__coverages is None:
            self.__coverages = np.fromiter(
                chain.from_iterable(
                    (coverage_for_gvcf_record(x) for x in self.records)
                ),
                # Coverages higher than 65535 will be incorrectly recorded
                # but this will save memory and these occurrences should be
                # very rare. Since coverage is never below 0 use an unsigned
                # integer.
                dtype=np.uint16
            )
        return self.__coverages

    @property
    def gq_qualities(self) -> np.ndarray:
        if self.__gq_qualities is None:
            self.__gq_qualities = np.fromiter(
                chain.from_iterable(
                    (gq_for_gvcf_record(x) for x in self.records)
                ),
                # GQ can never be higher than 99 and not lower than 0.
                # uint8 with values 0-255 is appropriate.
                dtype=np.uint8
            )
        return self.__gq_qualities

    @property
    def median_cov(self) -> float:
        return np.median(self.coverages)

    @property
    def mean_cov(self) -> float:
        return np.mean(self.coverages)

    @property
    def median_gq(self) -> float:
        return np.median(self.gq_qualities)

    @property
    def mean_gq(self) -> float:
        return qualmean(self.gq_qualities)

    @property
    def stats(self) -> dict:
        s = {
            "median_dp": self.median_cov,
            "mean_dp": self.mean_cov,
            "median_gq": self.median_gq,
            "mean_gq": self.mean_gq
        }
        coverage_boundaries = (10, 20, 30, 50, 100)
        coverage_at_least = fractions_at_least(self.coverages,
                                               coverage_boundaries)
        for perc, fraction in zip(coverage_boundaries, coverage_at_least):
            key = "perc_at_least_{0}_dp".format(perc)
            s[key] = fraction * 100

        gq_boundaries = (10, 20, 30, 50, 90)
        gq_at_least = fractions_at_least(self.gq_qualities,
                                         gq_boundaries)
        for perc, fraction in zip(gq_boundaries, gq_at_least):
            key = "perc_at_least_{0}_gq".format(perc)
            s[key] = fraction * 100
        return s


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


def vcf_records_to_coverage_array(vcf_records: Iterable[cyvcf2.Variant],
                                  maxlen: int = 15000
                                 ) -> np.ndarray:
    coverages = itertools.chain(coverage_for_gvcf_record(vcf_record, maxlen)
                                for vcf_record in vcf_records)
    return np.fromiter(coverages, dtype=np.uint16)


def region_coverages(reader: cyvcf2.VCF, regions: List[Region]) -> Optional[dict]:
    records = []
    for region in regions:
        it = reader(str(region))
        records += list(it)

    if len(records) == 0:
        return CovStats([]).stats

    return CovStats(records).stats
