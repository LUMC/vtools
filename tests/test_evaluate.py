import pytest

from vtools.evaluate import site_concordancy

from cyvcf2 import VCF


@pytest.fixture(scope='module')
def known_concordant():
    filename = 'tests/cases/gatk.vcf.gz'
    call = VCF(filename, gts012=True)
    positive = VCF(filename, gts012=True)
    d, disc = site_concordancy(call, positive, call_samples=['NA12878'],
                               positive_samples=['NA12878'],
                               min_gq=0, min_dp=0)
    return d


def test_total_sites(known_concordant):
    assert known_concordant['total_sites'] == 37


def test_sites_considered(known_concordant):
    assert known_concordant['sites_considered'] == 37


def test_alleles_considered(known_concordant):
    assert known_concordant['alleles_considered'] == 74


def test_alleles_het_concordant(known_concordant):
    assert known_concordant['alleles_het_concordant'] == 42


def test_alleles_hom_alt_concordant(known_concordant):
    assert known_concordant['alleles_hom_alt_concordant'] == 18


def test_alleles_hom_ref_concordant(known_concordant):
    assert known_concordant['alleles_hom_ref_concordant'] == 14


def test_alleles_concordant(known_concordant):
    assert known_concordant['alleles_concordant'] == 74


def test_alleles_discordant(known_concordant):
    assert known_concordant['alleles_discordant'] == 0


def test_alleles_no_call(known_concordant):
    assert known_concordant['alleles_no_call'] == 0


def test_alleles_low_qual(known_concordant):
    assert known_concordant['alleles_low_qual'] == 0


def test_alleles_low_depth(known_concordant):
    assert known_concordant['alleles_low_depth'] == 0


@pytest.fixture(scope='module')
def BLANK_NA12878():
    filename = 'tests/cases/gatk.vcf.gz'
    call = VCF(filename, gts012=True)
    positive = VCF(filename, gts012=True)
    d, disc = site_concordancy(call, positive, call_samples=['BLANK'],
                               positive_samples=['NA12878'],
                               min_gq=30, min_dp=20)
    return d


def test_low_qual_30(BLANK_NA12878):
    assert BLANK_NA12878['alleles_low_qual'] == 34


def test_low_depth_20(BLANK_NA12878):
    assert BLANK_NA12878['alleles_low_depth'] == 36


def test_no_call(BLANK_NA12878):
    assert BLANK_NA12878['alleles_no_call'] == 8


def test_truncated_call_file():
    """ When the call set is truncated, i.e. is missing variants which are
    present in the positive vcf file """

    filename = 'tests/cases/gatk.vcf.gz'
    truncated = 'tests/cases/gatk_truncated.vcf.gz'
    call = VCF(truncated, gts012=True)
    positive = VCF(filename, gts012=True)
    d, disc = site_concordancy(call, positive, call_samples=['NA12878'],
                               positive_samples=['BLANK'],
                               min_gq=30, min_dp=20)
    assert d['alleles_no_call'] == 12


def test_truncated_positive_file():
    """ When the positive set is truncated, i.e. the called vcf file contains
    variants which are absent from the positive vcf file """
    filename = 'tests/cases/gatk.vcf.gz'
    truncated = 'tests/cases/gatk_truncated.vcf.gz'
    call = VCF(filename, gts012=True)
    positive = VCF(truncated, gts012=True)
    d, disc = site_concordancy(call, positive, call_samples=['NA12878'],
                               positive_samples=['BLANK'],
                               min_gq=30, min_dp=20)
    assert d['alleles_no_call'] == 0


def test_phased_positive_file():
    """ Test error message when the positive vcf contains phased variants """
    filename = 'tests/cases/gatk.vcf.gz'
    phased = 'tests/cases/dummy_phased_blank.vcf.gz'
    call = VCF(filename, gts012=True)
    phased = VCF(phased, gts012=True)
    with pytest.raises(NotImplementedError,
                       match='Phased variants are not supported'):
        site_concordancy(call, phased, call_samples=['NA12878'],
                         positive_samples=['BLANK'])


def test_phased_call_file():
    """ Test error message when the call vcf contains phased variants """
    filename = 'tests/cases/gatk.vcf.gz'
    phased = 'tests/cases/dummy_phased_blank.vcf.gz'
    call = VCF(phased, gts012=True)
    positive = VCF(filename, gts012=True)
    with pytest.raises(NotImplementedError,
                       match='Phased variants are not supported'):
        site_concordancy(call, positive, call_samples=['BLANK'],
                         positive_samples=['NA12878'])


def test_decomposed_positive_file():
    """ Test error message when the positive vcf contains decomposed variants
    """
    filename = 'tests/cases/gatk.vcf.gz'
    decomposed = 'tests/cases/dummy_decomposed.vcf.gz'
    call = VCF(filename, gts012=True)
    positive = VCF(decomposed, gts012=True)
    with pytest.raises(NotImplementedError,
                       match='Decomposed variants are not supported'):
        site_concordancy(call, positive, call_samples=['BLANK'],
                         positive_samples=['NA12878'])


def test_decomposed_call_file():
    """ Test error message when the positive vcf contains decomposed variants
    """
    filename = 'tests/cases/gatk.vcf.gz'
    decomposed = 'tests/cases/dummy_decomposed.vcf.gz'
    call = VCF(decomposed, gts012=True)
    positive = VCF(filename, gts012=True)
    with pytest.raises(NotImplementedError,
                       match='Decomposed variants are not supported'):
        site_concordancy(call, positive, call_samples=['BLANK'],
                         positive_samples=['NA12878'])


def test_haploid_positive_file():
    """ Test error message when the positive vcf contains non-diploid variants
    """
    filename = 'tests/cases/gatk.vcf.gz'
    haploid = 'tests/cases/dummy_haploid.vcf.gz'
    call = VCF(filename, gts012=True)
    positive = VCF(haploid, gts012=True)
    with pytest.raises(NotImplementedError,
                       match='Non-diploid variants are not supported'):
        site_concordancy(call, positive, call_samples=['NA12878'],
                         positive_samples=['BLANK'])


def test_haploid_call_file():
    """ Test error message when the call vcf contains non-diploid variants
    """
    filename = 'tests/cases/gatk.vcf.gz'
    haploid = 'tests/cases/dummy_haploid.vcf.gz'
    call = VCF(haploid, gts012=True)
    positive = VCF(filename, gts012=True)
    with pytest.raises(NotImplementedError,
                       match='Non-diploid variants are not supported'):
        site_concordancy(call, positive, call_samples=['BLANK'],
                         positive_samples=['NA12878'])


@pytest.fixture(scope='module')
def partial_call_file():
    """ Test statistics when the call vcf contains partial variants """
    filename = 'tests/cases/gatk.vcf.gz'
    partial = 'tests/cases/gatk_partial_call.vcf.gz'
    call = VCF(partial, gts012=True)
    positive = VCF(filename, gts012=True)
    d, disc = site_concordancy(call, positive, call_samples=['BLANK'],
                               positive_samples=['BLANK'], min_dp=0, min_gq=0)
    return d


def test_partial_call_total_sites(partial_call_file):
    assert partial_call_file['total_sites'] == 37


def test_partial_call_alleles_hom_ref_concordant(partial_call_file):
    assert partial_call_file['alleles_hom_ref_concordant'] == 0


def test_partial_call_alleles_het_concordant(partial_call_file):
    assert partial_call_file['alleles_het_concordant'] == 0


def test_partial_call_alleles_hom_alt_concordant(partial_call_file):
    assert partial_call_file['alleles_hom_alt_concordant'] == 0


def test_partial_call_alleles_concordant(partial_call_file):
    assert partial_call_file['alleles_concordant'] == 6


def test_partial_call_alleles_discordant(partial_call_file):
    assert partial_call_file['alleles_discordant'] == 0


def test_partial_call_alleles_no_call(partial_call_file):
    assert partial_call_file['alleles_no_call'] == 68


@pytest.fixture(scope='module')
def partial_positive_file():
    """ Test statistics when the call vcf contains partial variants """
    filename = 'tests/cases/gatk.vcf.gz'
    partial = 'tests/cases/gatk_partial_call.vcf.gz'
    call = VCF(filename, gts012=True)
    positive = VCF(partial, gts012=True)
    d, disc = site_concordancy(call, positive, call_samples=['BLANK'],
                               positive_samples=['BLANK'], min_dp=0, min_gq=0)
    return d


def test_partial_positive_total_sites(partial_positive_file):
    assert partial_positive_file['total_sites'] == 6


def test_partial_positive_hom_ref_concordant(partial_positive_file):
    assert partial_positive_file['alleles_hom_ref_concordant'] == 0


def test_partial_positive_het_concordant(partial_positive_file):
    assert partial_positive_file['alleles_het_concordant'] == 0


def test_partial_positive_hom_alt_concordant(partial_positive_file):
    assert partial_positive_file['alleles_hom_alt_concordant'] == 0


def test_partial_positive_concordant(partial_positive_file):
    assert partial_positive_file['alleles_concordant'] == 6


def test_partial_positive_no_call(partial_positive_file):
    assert partial_positive_file['alleles_no_call'] == 6
