"""
vtools.evaluate
~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from typing import List, Dict
from types import SimpleNamespace

from cyvcf2 import VCF


def parse_variants(call: List[any], pos: List[any], results: Dict[str, int]):
    """ Parse the variants and add to results """

    call_variant = sorted(call[0:2])
    pos_variant = sorted(pos[0:2])

    # The types of concordant calls are counted separately
    if call_variant == pos_variant:
        # These variants are homozygous reference
        if call_variant == [0, 0]:
            results['alleles_hom_ref_concordant'] += 2
        # These variants are heterozygous, since they have different calls
        elif call_variant[0] != call_variant[1]:
            results['alleles_het_concordant'] += 2
        # If they are not homozygous reference, and not heterozygous, they must
        # be homozygous alternative, for whichever alt allele they have
        else:
            results['alleles_hom_alt_concordant'] += 2

    # Here we count all alleles independently, also for the concordant calls
    for i, j in zip(call_variant, pos_variant):
        if i == j:
            results['alleles_concordant'] += 1
        else:
            results['alleles_discordant'] += 1


def site_concordancy(call_vcf: VCF,
                     positive_vcf: VCF,
                     call_samples: List[str],
                     positive_samples: List[str],
                     min_gq: float = 30,
                     min_dp: float = 0) -> Dict[str, float]:
    """
    Calculate concordance between sites of two call sets,
    of which one contains known true positives.

    This only takes those sites in consideration which are shared between
    both call sets. This makes this evaluation useful in cases where
    the true positive set is known to contain regions not present at all
    in the base VCF. Such a situation occurs, for example, when comparing
    WES to a SNP array.

    Site identity is determined by the CHROM, POS, REF and ALT fields.
    This only works for bi-allelic variants.

    :param call_vcf: The VCF to evaluate. Must have gts012=True
    :param positive_vcf: The VCF file containing true positive sites.
        Must have gts012=True
    :param call_samples: List of samples in `call_vcf` to consider.
    :param positive_samples: List of samples in `positive_vcf` to consider.
    Must have same length than `call_samples`
    :param min_qg: Minimum quality of variants to consider.
    :param min_dp: Minimum depth of variants to consider.

    :raises: ValueError if sample lists are not of same size.

    :returns: Dictionary of shape {"hom_alt_concordant": 0,
                                   "het_concordant": 0, "sites_considered": 0,
                                   "total_sites": 0}
    """
    if len(positive_samples) != len(call_samples):
        raise ValueError("Lists of samples must have same size")

    pos_sample_idx = [positive_vcf.samples.index(x) for x in positive_samples]
    call_sample_idx = [call_vcf.samples.index(x) for x in call_samples]

    d = {
        "total_sites": 0,
        "sites_considered": 0,
        "alleles_considered": 0,
        "alleles_het_concordant": 0,
        "alleles_hom_alt_concordant": 0,
        "alleles_hom_ref_concordant": 0,
        "alleles_concordant": 0,
        "alleles_discordant": 0,
        "alleles_no_call": 0,
        "alleles_low_qual": 0,
        "alleles_low_depth": 0
    }

    # Keep track of the discordant sites
    discordant_count = 0
    discordant_records = list()

    # Keep track of the previous record, to make sure the positive vcf file
    # does not contain decomposed records
    prev_record = SimpleNamespace()
    prev_record.POS = -1
    prev_record.CHROM = -1

    for pos_record in positive_vcf:
        # Check to make sure the POS and CHROM are different from the previous
        # variant
        if (prev_record.POS == pos_record.POS and
                prev_record.CHROM == pos_record.CHROM):
            raise NotImplementedError('Decomposed variants are not supported')
        else:
            prev_record.POS = pos_record.POS
            prev_record.CHROM = pos_record.CHROM

        d['total_sites'] += 1
        query_str = "{0}:{1}-{2}".format(
            pos_record.CHROM,
            pos_record.start+1,
            pos_record.end
        )

        call_records = list(call_vcf(query_str))

        # If the variant is not present in the call vcf
        if len(call_records) == 0:
            d['alleles_no_call'] += 2
            continue

        # If there are multiple variants with the same CHROM and POS, the
        # variant in the call vcf file must have been decomposted
        if len(call_records) > 1:
            raise NotImplementedError('Decomposed variants are not supported')

        # Now we know there is only one call record
        call_record = call_records[0]

        d['sites_considered'] += 1
        d['alleles_considered'] += (2 * len(positive_samples))

        for p_s, c_s in zip(pos_sample_idx, call_sample_idx):
            p_gt = pos_record.gt_types[p_s]
            c_gt = call_record.gt_types[c_s]

            # Is there missing data in one of the vcf files?
            if p_gt == 3 or c_gt == 3:
                d['alleles_no_call'] += 2
                continue

            # Get the quality requirements for the call site
            c_gq = call_record.gt_quals[c_s]
            c_dp = call_record.gt_depths[c_s]

            # Did we fail any of the quality control checks?
            # We use this variable since we want to count all quality checks
            # independenly, and a single sample can fail multiple checks
            failed_quality_check = False

            # Is the QG of the call site below the threshold?
            if c_gq < min_gq:
                d['alleles_low_qual'] += 2
                failed_quality_check = True

            # Is the DP of the call site below the threshold?
            if c_dp < min_dp:
                d['alleles_low_depth'] += 2
                failed_quality_check = True

            # Go to the next variant if any of the quality checks failed
            if failed_quality_check:
                continue

            # Get the genotypes
            pos = pos_record.genotypes[p_s]
            cal = call_record.genotypes[c_s]

            # If the genotypes are phased
            if pos[2] or cal[2]:
                raise NotImplementedError('Phased variants are not supported')

            # Parse the genotypes and add the results into d
            parse_variants(pos, cal, d)

            # The current variant is discordant
            if d['alleles_discordant'] > discordant_count:
                discordant_count = d['alleles_discordant']
                discordant_records.append(call_record)

    return d, discordant_records
