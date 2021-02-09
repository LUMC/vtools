"""
vtools.stats
~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden Unversity Medical Center
:license: MIT
"""

import json
from collections import Counter

import cyvcf2  # type: ignore

from .optimized import Sample  # type: ignore


def gen_chrom_counter(vcf_reader):
    """Generate chromosome counter from VCF reader"""
    return Counter({n: 0 for n in vcf_reader.seqnames})


class Stats(object):
    def __init__(self, vcf_path):
        self.path = vcf_path
        self.vcf = cyvcf2.VCF(vcf_path, gts012=True)

        self.samples = [Sample(x, i) for i, x in enumerate(self.vcf.samples)]
        self.chrom_counter = gen_chrom_counter(self.vcf)

        self.__calculated = False

    def calculate(self):
        for record in self.vcf:
            for s in self.samples:
                s.add_variant(record)
            self.chrom_counter[record.CHROM] += 1
        self.__calculated = True

    @property
    def total_variants(self):
        return sum(self.chrom_counter.values())

    @property
    def as_dict(self):
        if not self.__calculated:
            self.calculate()
        return {
            "vcf_path": self.path,
            "total_variants": self.total_variants,
            "samples": [s.as_dict for s in self.samples],
            "per_chromosome_variants": {k: v for k, v in
                                        self.chrom_counter.items()}
        }

    @property
    def as_json(self):
        return json.dumps(self.as_dict, sort_keys=True, indent=4)
