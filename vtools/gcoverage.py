"""
vtools.gcoverage
~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

import cyvcf2
import numpy as np

from collections import namedtuple
from itertools import chain

from typing import List, Optional, Tuple

from .optimized import amount_atleast


Region = namedtuple("Region", ["chr", "start", "end"])


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
                dtype=int
            )
        return self.__coverages

    @property
    def gq_qualities(self) -> np.ndarray:
        if self.__gq_qualities is None:
            self.__gq_qualities = np.fromiter(
                chain.from_iterable(
                    (gq_for_gvcf_record(x) for x in self.records)
                ),
                dtype=int
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

    def percent_atleast_dp(self, atleast) -> Optional[float]:
        if len(self.coverages) == 0:
            return None
        k = amount_atleast(self.coverages, atleast)
        return (k/len(self.coverages))*100

    def percent_atleast_gq(self, atleast) -> Optional[float]:
        if len(self.gq_qualities) == 0:
            return None
        k = amount_atleast(self.gq_qualities, atleast)
        return (k/len(self.gq_qualities))*100

    @property
    def stats(self) -> dict:
        s = {
            "median_dp": self.median_cov,
            "mean_dp": self.mean_cov,
            "median_gq": self.median_gq,
            "mean_gq": self.mean_gq
        }
        for perc in (10, 20, 30, 50, 100):
            key = "perc_at_least_{0}_dp".format(perc)
            s[key] = self.percent_atleast_dp(perc)

        for perc in (10, 20, 30, 50, 90):
            key = "perc_at_least_{0}_gq".format(perc)
            s[key] = self.percent_atleast_gq(perc)
        return s


class RefRecord(object):
    def __init__(self, line: str):
        self.line = line  # type: str
        self.gene = None  # type: Optional[str]
        self.transcript = None  # type: Optional[str]
        self.contig = None  # type: Optional[str]
        self.start = None  # type: Optional[int]
        self.end = None  # type: Optional[int]
        self.cds_start = None  # type: Optional[int]
        self.cds_end = None  # type: Optional[int]
        self.exon_starts = []  # type: List[int]
        self.exon_ends = []  # type: List[int]
        self.forward = True

        self.parse()

    def parse(self):
        contents = self.line.strip().split("\t")
        if len(contents) < 11:
            raise ValueError("refFlat line must have at least 11 fields")
        self.gene = contents[0]
        self.transcript = contents[1]
        self.contig = contents[2]
        if "-" in contents[3].strip():
            self.forward = False
        self.start = int(contents[4])
        self.end = int(contents[5])
        self.cds_start = int(contents[6])
        self.cds_end = int(contents[7])
        self.exon_starts = [int(x) for x in contents[9].split(",")[:-1]]
        self.exon_ends = [int(x) for x in contents[10].split(",")[:-1]]

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


def region_coverages(reader: cyvcf2.VCF, regions: List[Region]) -> Optional[dict]:
    records = []
    for region in regions:
        reg_str = "{0}:{1}-{2}".format(region.chr, region.start, region.end)
        it = reader(reg_str)
        records += list(it)

    if len(records) == 0:
        return CovStats([]).stats

    return CovStats(records).stats
