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


def parse_variants(ref: str, call: List[str], pos: List[str],
                   results: Dict[str, int]):
    """ Parse the variants and add to results """

    call_variant = set(call)
    pos_variant = set(pos)

    # The types of concordant calls are counted separately
    if call_variant == pos_variant:
        # These variants are homozygous reference
        if len(call_variant) == 1 and next(iter(call_variant)) == ref:
            results['alleles_hom_ref_concordant'] += 2
        # These variants are heterozygous, since they have different calls
        elif len(call_variant) > 1:
            results['alleles_het_concordant'] += 2
        # If they are not homozygous reference, and not heterozygous, they must
        # be homozygous alternative, for whichever alt allele they have
        else:
            results['alleles_hom_alt_concordant'] += 2

    # Here we count all alleles independently, also for the concordant calls
    for allele in call:
        if allele == '.':
            results['alleles_no_call'] += 1
        elif allele in pos:
            results['alleles_concordant'] += 1
            # We cannot match an A/A call with A/G, so we need to remove the
            # A call from the positive set once we have 'used' it
            pos.remove(allele)
        else:
            results['alleles_discordant'] += 1


def RGQ_header_defined(vcf):
    """ Determine whether the RGQ annotation is defined in the FORMAT header of
    the vcf.
    """
    try:
        vcf.get_header_type('RGQ')
    except KeyError:
        return False
    else:
        return True


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

    # Determine if the 'RQC' annotation is present in the FORMAT header
    # This is needed for avoid getting an KeyError when we want to see if the
    # RQC annotation is set for a given variant
    RGQ_present = RGQ_header_defined(call_vcf)

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

            # Is one of the sites no call. This is only for fully no call
            # sites. Partial no call site are handled in parse_variants
            if p_gt == 3 or c_gt == 3:
                d['alleles_no_call'] += 2
                continue

            # Get the quality requirements for the call site
            c_dp = call_record.gt_depths[c_s]
            # If the 'RGQ' is defined, which is the case for homref calls from
            # GATK, use this value instead of 'GQ'
            # See
            # https://gatkforums.broadinstitute.org/gatk/discussion/9907/genotypegvcfs-no-records-in-vcf  # noqa
            if RGQ_present and call_record.format('RGQ') is not None:
                c_gq = call_record.format('RGQ')[0][0]
            else:
                c_gq = call_record.gt_quals[c_s]

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

            # If the genotypes are not diploid
            if len(pos) != 3 or len(cal) != 3:
                raise NotImplementedError('Non-diploid variants are not '
                                          'supported')
            # If the genotypes are phased
            if pos[2] or cal[2]:
                raise NotImplementedError('Phased variants are not supported')

            # Get the genotyped bases. There is a deprecationWarning when using
            # gt_bases directly so we go via the alleles and the gt call
            pos_alleles = [pos_record.REF] + pos_record.ALT
            cal_alleles = [call_record.REF] + call_record.ALT

            # The third item in pos is a boolean indicating phasing, so we
            # don't use that
            # No call is indicated by -1 in pyvcf2, so here we convert negative
            # values to '.', which is what is used in the vcf to indicate no
            # call
            pos_gt = [pos_alleles[x] if x >= 0 else '.' for x in pos[:2]]
            cal_gt = [cal_alleles[x] if x >= 0 else '.' for x in cal[:2]]

            # Parse the genotypes and add the results into d
            # We also need to know the reference call to determine if variants
            # are hom ref. For this, we take the positive vcf REF call to be
            # the truth
            parse_variants(pos_record.REF, cal_gt, pos_gt, d)

            # The current variant is discordant
            if d['alleles_discordant'] > discordant_count:
                discordant_count = d['alleles_discordant']
                discordant_records.append(call_record)

    return d, discordant_records
