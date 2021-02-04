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

from typing import List, Optional, Tuple, NamedTuple

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


def region_coverages(reader: cyvcf2.VCF, regions: List[Region]) -> Optional[dict]:
    records = []
    for region in regions:
        reg_str = "{0}:{1}-{2}".format(region.chr, region.start, region.end)
        it = reader(reg_str)
        records += list(it)

    if len(records) == 0:
        return CovStats([]).stats

    return CovStats(records).stats
