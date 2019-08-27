import pytest

from vtools.evaluate import site_concordancy

from cyvcf2 import VCF


@pytest.fixture
def known_concordant():
    known = 'tests/cases/gatk.vcf.gz'
    d, disc = site_concordancy(VCF(known, gts012=True), VCF(known,
                                                            gts012=True),
                               ['NA12878'], ['NA12878'], min_gq=0,
                               min_dp=0)
    return d


@pytest.fixture
def blank_NA12878():
    filename = 'tests/cases/gatk.vcf.gz'
    d, disc = site_concordancy(VCF(filename, gts012=True), VCF(filename,
                                                               gts012=True),
                               ['NA12878'], ['BLANK'], min_gq=0, min_dp=0)
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


def test_alleles_no_call(blank_NA12878):
    assert blank_NA12878['alleles_no_call'] == 8


def test_alleles_low_qual(known_concordant):
    assert known_concordant['alleles_low_qual'] == 0


def test_alleles_low_depth(known_concordant):
    assert known_concordant['alleles_low_depth'] == 0
