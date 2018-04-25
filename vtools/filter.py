"""
vtools.filter
~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from cyvcf2 import VCF, Writer

import json
import enum

from .optimized import get_af

class FilterClass(enum.Enum):
    NON_CANONICAL = {
        "ID": "NON_CANONICAL",
        "Description": "Non-canonical chromosome"
    }
    INDEX_UNCALLED = {
        "ID": "INDEX_UNCALLED",
        "Description": "Index uncalled or homozygous reference"
    }
    TOO_HIGH_GONL_AF = {
        "ID": "TOO_HIGH_GONL_AF",
        "Description": "Too high GoNL allele frequency"
    }
    TOO_HIGH_GNOMAD_AF = {
        "ID": "TOO_HIGH_GNOMAD_AF",
        "Description": "Too high GnomAD allele frequency"
    }
    LOW_GQ = {
        "ID": "LOW_GQ",
        "Description": "Too low GQ on index"
    }
    DELETED_ALLELE = {
        "ID": "DELETED_ALLELE",
        "Description": "The only ALT allele is a deleted allele."
    }


class FilterParams(object):

    def __init__(self, params_path):
        self.params_path = params_path

        with open(self.params_path) as handle:
            self.params_dict = json.load(handle)

        self.__gonl_vcf = None
        self.__onekg_vcf = None
        self.__gnomad_vcf = None

    @property
    def min_gq(self):
        return self.params_dict['gq_pass']

    @property
    def max_gonl_af(self):
        return self.params_dict['low_gonl_af']

    @property
    def max_gnomad_af(self):
        return self.params_dict['low_gnomad_af']

    @property
    def index_called(self):
        return self.params_dict['index_called']

    @property
    def canonical_chromosomes(self):
        return self.params_dict['canonical_chromosomes']

    @property
    def gonl_vcf(self):
        if self.__gonl_vcf is None:
            self.__gonl_vcf = VCF(self.params_dict['gonl_vcf'], gts012=True)
        return self.__gonl_vcf

    @property
    def gonl_af(self):
        return self.params_dict["gonl_af"]

    @property
    def gnomad_vcf(self):
        if self.__gnomad_vcf is None:
            self.__gnomad_vcf = VCF(self.params_dict['gnomad_vcf'], gts012=True)
        return self.__gnomad_vcf

    @property
    def gnomad_af(self):
        return self.params_dict['gnomad_af']


class Filterer(object):
    """
    Iterator returning tuples of:
    (cyvcf2_record, filter)

    if the filter item is None, record is unfiltered

    We assume gts012 to be set to True in the cyvcf2 readers
    """

    def __init__(self, vcf_it, filter_params, index, immediate_return=True):
        self.vcf_it = vcf_it
        self.filters = filter_params
        self.index = index
        self.immediate_return = immediate_return

        self.canonical_chroms = {"M", "X", "Y"}.union(set(map(str, range(0, 23))))

    def __next__(self):
        record = next(self.vcf_it)
        filters = []

        # filter out deleted alleles regardless
        if set(record.ALT) == {"*"}:
            return record, [FilterClass.DELETED_ALLELE]

        if self.filters.canonical_chromosomes:
            chrom = str(record.CHROM)
            if chrom.startswith("chr"):
                chrom = chrom.split("chr")[-1]
            if chrom not in self.canonical_chroms:
                if self.immediate_return:
                    return record, [FilterClass.NON_CANONICAL]
                filters.append(FilterClass.NON_CANONICAL)

        if self.filters.index_called:
            gt = record.gt_types[self.index]
            if gt == 0 or gt == 3:
                if self.immediate_return:
                    return record, [FilterClass.INDEX_UNCALLED]
                filters.append(FilterClass.INDEX_UNCALLED)

        if self.filters.min_gq is not None:
            gq = record.gt_quals[self.index]
            if gq < self.filters.min_gq:
                if self.immediate_return:
                    return record, [FilterClass.LOW_GQ]
                filters.append(FilterClass.LOW_GQ)

        if self.filters.max_gonl_af is not None:
            gonl_af = get_af(record, self.filters.gonl_vcf,
                                  self.filters.gonl_af)
            if gonl_af > self.filters.max_gonl_af:
                if self.immediate_return:
                    return record, [FilterClass.TOO_HIGH_GONL_AF]
                filters.append(FilterClass.TOO_HIGH_GONL_AF)

        if self.filters.max_gnomad_af is not None:
            gnomad_af = get_af(record, self.filters.gnomad_vcf,
                                  self.filters.gnomad_af)
            if gnomad_af > self.filters.max_gnomad_af:
                if self.immediate_return:
                    return record, [FilterClass.TOO_HIGH_GNOMAD_AF]
                filters.append(FilterClass.TOO_HIGH_GNOMAD_AF)

        return record, filters

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self
