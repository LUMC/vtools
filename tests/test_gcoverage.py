# MIT License
#
# Copyright (c) 2018, 2020 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import math
import sys
from pathlib import Path

from cyvcf2 import VCF  # type: ignore

import numpy as np

import pytest

from vtools import gcoverage
from vtools.gcoverage import CovStats, RefRecord,  Region, \
    feature_to_coverage_and_quality_lists, file_to_refflat_records, \
    qualmean, refflat_and_gvcfs_to_tsv, \
    region_and_vcf_to_coverage_and_quality_lists


DATA_DIR = Path(__file__).parent / "gcoverage_data"
TEST_GVCF = DATA_DIR / "test.g.vcf.gz"
TEST_REFFLAT = DATA_DIR / "test.refflat"


def test_qualmean():
    qualities = np.array([5, 15, 25, 35, 45, 75, 95])
    assert math.isclose(qualmean(qualities), 12.9934485294235)


def test_covstats_from_coverages_and_gq_qualities():
    # Normally coverages and qualities should have the same length but
    # these values are taken for easy testing.
    # The numbers in coverages are  randomly generated but are are sorted for
    # convenience: list(sorted([random.randint(0, 100) for _ in range(100)]))
    coverages = [
        1, 2, 7, 8, 8, 9, 9, 9, 9, 11, 11, 14, 17, 19, 22, 23, 23, 24, 24, 25,
        26, 27, 27, 27, 27, 29, 29, 32, 34, 34, 35, 36, 37, 39, 39, 41, 41, 44,
        44, 44, 46, 47, 47, 48, 48, 50, 50, 51, 53, 54, 54, 56, 56, 56, 57, 58,
        59, 59, 63, 63, 65, 66, 67, 67, 67, 70, 70, 71, 73, 73, 76, 76, 77, 77,
        79, 79, 80, 81, 82, 84, 84, 85, 86, 86, 86, 91, 91, 93, 93, 93, 93, 94,
        94, 95, 95, 95, 96, 96, 97, 97]
    qualities = [5, 15, 25, 35, 45, 75, 95]
    covstats = CovStats.from_coverages_and_gq_qualities(coverages, qualities)
    assert math.isclose(covstats.mean_gq, 12.9934485294235)
    assert covstats.median_gq == 35
    assert covstats.mean_dp == sum(coverages) / len(coverages)
    assert covstats.median_dp == 54
    assert covstats.perc_at_least_10_dp == 91 / 100 * 100
    assert covstats.perc_at_least_20_dp == 86 / 100 * 100
    assert covstats.perc_at_least_30_dp == 73 / 100 * 100
    assert covstats.perc_at_least_50_dp == 55 / 100 * 100
    assert covstats.perc_at_least_100_dp == 0 / 100 * 100
    assert covstats.perc_at_least_10_gq == 6 / 7 * 100
    assert covstats.perc_at_least_20_gq == 5 / 7 * 100
    assert covstats.perc_at_least_30_gq == 4 / 7 * 100
    assert covstats.perc_at_least_50_gq == 2 / 7 * 100
    assert covstats.perc_at_least_90_gq == 1 / 7 * 100


def test_refrecord():
    record_line = ("GENE\tTRANSCRIPT\tcontig\t-\t100\t1000\t300\t600\t5\t"
                   "100,250,400,550,800,\t200,350,500,650,900,")
    refflat_record = RefRecord.from_line(record_line)
    assert refflat_record.gene == "GENE"
    assert refflat_record.transcript == "TRANSCRIPT"
    assert refflat_record.contig == "contig"
    assert refflat_record.forward is False
    assert refflat_record.start == 100
    assert refflat_record.end == 1000
    assert refflat_record.cds_start == 300
    assert refflat_record.cds_end == 600
    assert refflat_record.exon_starts == [100, 250, 400, 550, 800]
    assert refflat_record.exon_ends == [200, 350, 500, 650, 900]
    assert refflat_record.exons == [
        Region("contig", 100, 200),
        Region("contig", 250, 350),
        Region("contig", 400, 500),
        Region("contig", 550, 650),
        Region("contig", 800, 900)
    ]
    # Only (partial) exons inside the cds_start - cds_end boundary should end
    # up here.
    assert refflat_record.cds_exons == [
        Region("contig", 300, 350),
        Region("contig", 400, 500),
        Region("contig", 550, 600)
    ]


def test_refrecord_cds_exons():
    refrecord_line = ("GENE1\tTR0001\tchr1\t+\t11\t20\t13\t17\t3\t"
                      "11,14,17,\t12,16,19,")
    refrecord = RefRecord.from_line(refrecord_line)
    assert len(refrecord.cds_exons) == 2
    assert refrecord.cds_exons == [Region("chr1", 14, 16),
                                   Region("chr1", 17, 17)]


def test_refrecord_invalid_strand():
    record_line = ("GENE\tTRANSCRIPT\tcontig\tPLUS\t100\t1000\t300\t600\t5\t"
                   "100,250,400,550,800,\t200,350,500,650,900,")
    with pytest.raises(ValueError) as error:
        RefRecord.from_line(record_line)
    error.match("Invalid strand")


def test_refrecord_forward():
    record_line = ("GENE\tTRANSCRIPT\tcontig\t+\t100\t1000\t300\t600\t5\t"
                   "100,250,400,550,800,\t200,350,500,650,900,")
    refflat_record = RefRecord.from_line(record_line)
    assert refflat_record.forward is True


def test_refrecord_too_many_columns():
    record_line = ("GENE\tTRANSCRIPT\tcontig\t+\t100\t1000\t300\t600\t5\t"
                   "100,250,400,550,800,\t200,350,500,650,900,\tnonsense")
    with pytest.raises(ValueError) as error:
        RefRecord.from_line(record_line)
    error.match("exactly 11")


def test_refrecord_too_little_columns():
    record_line = ("GENE\t+\t100\t1000\t300\t600\t5\t"
                   "100,250,400,550,800,\t200,350,500,650,900,\tnonsense")
    with pytest.raises(ValueError) as error:
        RefRecord.from_line(record_line)
    error.match("exactly 11")


def test_region_string():
    assert str(Region("chr1", 1, 100)) == "chr1:1-100"


def test_region_length():
    assert len(Region("chr1", 1, 100)) == 100


def test_file_ro_refflat_records():
    test_refflat = DATA_DIR / "10genes.refflat"
    record_list = list(file_to_refflat_records(test_refflat))
    transcripts = [record.transcript for record in record_list]
    assert transcripts == [
        "NR_026971", "NM_000014", "NM_001173466", "NM_015665", "NM_001271885",
        "NM_001271886", "NM_024666", "NM_020745", "NM_001605", "NM_005763"]


def test_region_and_vcf_to_coverage_and_quality_lists():
    test_vcf = VCF(str(TEST_GVCF))
    region = Region("chr1", 20, 30)
    coverages, qualities = region_and_vcf_to_coverage_and_quality_lists(
        region, test_vcf)
    assert len(coverages) == len(region)
    assert len(coverages) == len(qualities)
    assert coverages == [7, 8, 9, 9, 9, 10, 10, 10, 10, 10, 10]
    assert qualities == [21.0, 24.0, 27.0, 27.0, 27.0, 30.0, 30.0, 30.0, 30.0,
                         30.0, 30.0]


def test_feature_to_coverage_and_quality_arrays():
    test_vcf = VCF(str(TEST_GVCF))
    regions = [Region("chr1", 1, 10), Region("chr1", 21, 30)]
    coverages, qualities = feature_to_coverage_and_quality_lists(
        regions, [test_vcf])
    assert len(coverages) == len(qualities)
    assert len(coverages) == sum(len(region) for region in regions)
    assert coverages == [2, 3, 4, 4, 4, 5, 5, 6, 6, 6, 8, 9, 9, 9, 10, 10, 10,
                         10, 10, 10]
    assert qualities == [6.0, 9.0, 12.0, 12.0, 12.0, 15.0, 15.0, 18.0, 0.0,
                         18.0, 24.0, 27.0, 27.0, 27.0, 30.0, 30.0, 30.0, 30.0,
                         30.0, 30.0]


def test_refflat_and_gvcfs_to_tsv_per_exon_compact():
    lines = refflat_and_gvcfs_to_tsv(TEST_REFFLAT, [TEST_GVCF],
                                     region_of_interest="exon",
                                     compact_header=True)
    result = list(lines)
    assert result[0] == ("gene\ttranscript\texon\t" + CovStats.header(True))
    assert result[1] == ("GENE1\tTR0001\t1\t6.00\t18.00\t6.00\t18.00\t"
                         "0.00\t0.00\t0.00\t0.00\t0.00\t"
                         "100.00\t0.00\t0.00\t0.00\t0.00")


def test_refflat_and_gvcfs_to_tsv_per_transcript_verbose():
    lines = refflat_and_gvcfs_to_tsv(TEST_REFFLAT, [TEST_GVCF],
                                     region_of_interest="transcript",
                                     compact_header=False)
    result = list(lines)
    assert result[0] == ("gene\ttranscript\t" + CovStats.header(False))
    assert result[3] == ("GENE3\tTR0006\t32.40\t49.23\t32.00\t79.00\t"
                         "100.00\t100.00\t100.00\t0.00\t0.00\t"
                         "100.00\t100.00\t100.00\t80.00\t20.00")


def test_refflat_and_gvcfs_to_tsv_per_transcript_cds_exons():
    lines = refflat_and_gvcfs_to_tsv(TEST_REFFLAT, [TEST_GVCF],
                                     region_of_interest="transcript_cds_exons",
                                     compact_header=False)
    result = list(lines)
    assert result[0] == ("gene\ttranscript\t" + CovStats.header(False))
    assert result[1].startswith("GENE1\tTR0001")
    assert result[2].startswith("GENE3\tTR0006")
    # GENE2 Should be omitted as it has now coding exons
    assert len(result) == 3
    assert result[1] == ("GENE1\tTR0001\t6.00\t5.82\t6.00\t18.00\t"
                         "0.00\t0.00\t0.00\t0.00\t0.00\t"
                         "75.00\t0.00\t0.00\t0.00\t0.00")


MAIN_TESTS = [
    (["-R", str(TEST_REFFLAT), str(TEST_GVCF)],
     "GENE1\tTR0001\t1\t6.00\t18.00\t6.00\t18.00\t0.00\t0.00\t0.00\t0.00\t"
     "0.00\t100.00\t0.00\t0.00\t0.00\t0.00"
     ),
    (["-R", str(TEST_REFFLAT), str(TEST_GVCF), "--per-transcript"],
     "GENE3\tTR0006\t32.40\t49.23\t32.00\t79.00\t100.00\t100.00\t100.00\t"
     "0.00\t0.00\t100.00\t100.00\t100.00\t80.00\t20.00"
     ),
    (["-R", str(TEST_REFFLAT), str(TEST_GVCF),
      "--per-transcript-cds-exons"],
     "GENE1\tTR0001\t6.00\t5.82\t6.00\t18.00\t0.00\t0.00\t0.00\t0.00\t"
     "0.00\t75.00\t0.00\t0.00\t0.00\t0.00"
     )
]


@pytest.mark.parametrize(["args", "in_result"], MAIN_TESTS)
def test_main(args, in_result, capsys):
    sys.argv = [""] + args
    gcoverage.main()
    result = capsys.readouterr()
    assert in_result in result.out
