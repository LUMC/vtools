"""
vtools.optimized
~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from collections import Counter
import numpy
cimport numpy as np

cpdef int amount_atleast(np.int64_t[::1] values, int atleast):
    """
    Return amount of values at least `atleast`
    :param values: Iterable of int
    :param atleast: int
    :return: int
    """
    cdef int passed = 0
    cdef int val
    cdef size_t i
    for i in range(values.shape[0]):
        val = values[i]
        if val >= atleast:
            passed += 1
    return passed


cdef class Sample(object):
    cdef unicode name
    cdef size_t idx
    cdef int transversions
    cdef int transitions
    cdef int hom_ref
    cdef int het
    cdef int hom_alt
    cdef int deletions
    cdef int insertions
    cdef int snps
    cdef int gq_counter[100]
    def __init__(self, name, size_t idx):
        self.name = name
        self.idx = idx
        self.transversions = 0
        self.transitions = 0
        self.hom_ref = 0
        self.het = 0
        self.hom_alt = 0
        self.deletions = 0
        self.insertions = 0
        self.snps = 0

    @property
    def ti_tv(self):
        cdef float titv
        if self.transversions > 0:
            titv = float(self.transitions)/self.transversions
        else:
            titv = numpy.nan
        return titv

    cpdef void add_variant(self, var):
        """assuming gts012=True"""
        cdef int typ = var.gt_types[self.idx]
        if typ == 3:
            return

        cdef int gq = var.gt_quals[self.idx]
        self.gq_counter[gq] += 1

        if typ == 0:
            self.hom_ref += 1
            return

        if typ == 1:
            self.het += 1
        if typ == 2:
            self.hom_alt += 1

        if var.is_snp and var.is_transition:
            self.transitions += 1
            self.snps += 1
        elif var.is_snp:
            self.transversions += 1
            self.snps += 1
        elif var.is_indel and var.is_deletion:
            self.deletions += 1
        elif var.is_indel:
            self.insertions += 1
        return

    @property
    def gq_distr(self):
        return [self.gq_counter[x] for x in range(100)]

    @property
    def total_variants(self):
        return self.hom_alt + self.het

    @property
    def as_dict(self):
        return {
            "name": self.name,
            "total_variants": self.total_variants,
            "variant_types": {
                "snps": self.snps,
                "deletions": self.deletions,
                "insertions": self.insertions,
            },
            "genotypes": {
                "hom_ref": self.hom_ref,
                "het": self.het,
                "hom_alt": self.hom_alt
            },
            "transitions": self.transitions,
            "transversions": self.transversions,
            "ti_tv_ratio": self.ti_tv,
            "gq_distribution": self.gq_distr
        }

