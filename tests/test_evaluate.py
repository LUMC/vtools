import pytest

from vtools.evaluate import site_concordancy

from cyvcf2 import VCF


@pytest.fixture
def known_concordant():
    known = '/home/rrvandenberg/devel/vtools/tests/cases/gatk.vcf.gz'
    d, disc = site_concordancy(VCF(known, gts012=True), VCF(known,
                                                            gts012=True),
                               ['NA12878'], ['NA12878'], min_gq=0,
                               min_dp=0)
    return d


def test_total_sites(known_concordant):
    assert known_concordant['total_sites'] == 37


def test_sites_considered(known_concordant):
    assert known_concordant['sites_considered'] == 37


def test_alleles_considered(known_concordant):
    assert known_concordant['alleles_considered'] == 74


def test_alleles_het_concordant(known_concordant):
    assert known_concordant['alleles_het_concordant'] == 42
