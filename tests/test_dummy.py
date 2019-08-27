import pytest

from vtools.evaluate import site_concordancy

from cyvcf2 import VCF


@pytest.fixture
def known():
    bla =VCF('tests/cases/known1.vcf.gz', gts012=True)
    return bla

def test_evaluate(known):
    d, disc = site_concordancy(known, known, ["known1"], ["known1"])
    assert d['total_sites'] == 3
    assert not disc