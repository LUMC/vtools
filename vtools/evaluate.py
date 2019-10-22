"""
vtools.evaluate
~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from typing import List, Dict

from cyvcf2 import VCF


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
    discordant_count = 0
    discordant_records = list()

    for pos_record in positive_vcf:
        d['total_sites'] += 1
        query_str = "{0}:{1}-{2}".format(
            pos_record.CHROM,
            pos_record.start,
            pos_record.end
        )

        it = call_vcf(query_str)
        same = []
        # If the vcf file has been decomposed, there can be multiple sites with
        # the same CHROM and POS, which is why we have to iterate over all the
        # sites that are returned
        for it_record in it:
            # We only want to consider sites that have the same CHROM, POS, REF
            # and ALT as the positive vcf, since otherwise we might consider
            # A/T and A/G as identical since they are both heterozygous
            if (it_record.CHROM == pos_record.CHROM
                    and it_record.POS == pos_record.POS
                    and it_record.REF == pos_record.REF
                    and it_record.ALT == pos_record.ALT):
                same.append(it_record)

        # If the variant is not present in the call vcf
        if len(same) == 0:
            d['alleles_no_call'] += 2

        if len(same) != 1:
            continue

        call_record = same[0]
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

            # If the site passed the quality requirements
            if p_gt == 0 and c_gt == 0:
                # both homozygous reference
                d['alleles_concordant'] += 2
                d['alleles_hom_ref_concordant'] += 2
            elif p_gt == 0 and c_gt == 1:
                # heterozygous when homozygous reference
                d['alleles_concordant'] += 1
                d['alleles_discordant'] += 1
            elif p_gt == 0 and c_gt == 2:
                # homozygous alt when homozygous ref
                d['alleles_discordant'] += 2
            elif p_gt == 1 and c_gt == 0:
                # homozygous ref when heterozygous
                d['alleles_concordant'] += 1
                d['alleles_discordant'] += 1
            elif p_gt == 1 and c_gt == 1:
                # both heterozygous
                d['alleles_concordant'] += 2
                d['alleles_het_concordant'] += 2
            elif p_gt == 1 and c_gt == 2:
                # homozygous alt when heterozygous
                d['alleles_concordant'] += 1
                d['alleles_discordant'] += 1
            elif p_gt == 2 and c_gt == 0:
                # homozygous ref when homozygous alt
                d['alleles_discordant'] += 2
            elif p_gt == 2 and c_gt == 1:
                # heterozygous when homozygous alt
                d['alleles_concordant'] += 1
                d['alleles_discordant'] += 1
            elif p_gt == 2 and c_gt == 2:
                # both homozygous alt
                d['alleles_concordant'] += 2
                d['alleles_hom_alt_concordant'] += 2
            else:
                # This should not be reached
                raise RuntimeError('Unhandled genotype call')

            # The current variant is discordant
            if d['alleles_discordant'] > discordant_count:
                discordant_count = d['alleles_discordant']
                discordant_records.append(call_record)

    return d, discordant_records
