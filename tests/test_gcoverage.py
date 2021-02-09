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
from pathlib import Path

import numpy as np

import pytest

from vtools.gcoverage import CovStats, RefRecord,  Region, \
    file_to_refflat_records, qualmean


def test_qualmean():
    qualities = np.array([5, 15, 25, 35, 45, 75, 95])
    assert math.isclose(qualmean(qualities), 12.9934485294235)


def test_covstats_from_coverages_and_gq_qualities():
    # Normally coverages and qualities should have the same length but
    # these values are taken for easy testing.
    # The numbers in coverages are  randomly generated but are are sorted for
    # convenience: list(sorted([random.randint(0, 100) for _ in range(100)]))
    coverages = np.array([
        1, 2, 7, 8, 8, 9, 9, 9, 9, 11, 11, 14, 17, 19, 22, 23, 23, 24, 24, 25,
        26, 27, 27, 27, 27, 29, 29, 32, 34, 34, 35, 36, 37, 39, 39, 41, 41, 44,
        44, 44, 46, 47, 47, 48, 48, 50, 50, 51, 53, 54, 54, 56, 56, 56, 57, 58,
        59, 59, 63, 63, 65, 66, 67, 67, 67, 70, 70, 71, 73, 73, 76, 76, 77, 77,
        79, 79, 80, 81, 82, 84, 84, 85, 86, 86, 86, 91, 91, 93, 93, 93, 93, 94,
        94, 95, 95, 95, 96, 96, 97, 97])
    qualities = np.array([5, 15, 25, 35, 45, 75, 95])
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


def test_refrecord_forward():
    record_line = ("GENE\tTRANSCRIPT\tcontig\t+\t100\t1000\t300\t600\t5\t"
                   "100,250,400,550,800,\t200,350,500,650,900,")
    refflat_record = RefRecord.from_line(record_line)
    assert refflat_record.forward is True


def test_refrecord_too_many_columns():
    record_line = ("GENE\tTRANSCRIPT\tcontig\t+\t100\t1000\t300\t600\t5\t"
                   "100,250,400,550,800,\t200,350,500,650,900,\tnonsense")
    with pytest.raises(ValueError):
        RefRecord.from_line(record_line)


def test_refrecord_too_little_columns():
    record_line = ("GENE\t+\t100\t1000\t300\t600\t5\t"
                   "100,250,400,550,800,\t200,350,500,650,900,\tnonsense")
    with pytest.raises(ValueError):
        RefRecord.from_line(record_line)


def test_region_string():
    assert str(Region("chr1", 1, 100)) == "chr1:1-100"


def test_file_ro_refflat_records():
    test_refflat = Path(__file__).parent / "gcoverage_data" / "10genes.refflat"
    record_list = list(file_to_refflat_records(test_refflat))
    transcripts = [record.transcript for record in record_list]
    assert transcripts == [
        "NR_026971", "NM_000014", "NM_001173466", "NM_015665", "NM_001271885",
        "NM_001271886", "NM_024666", "NM_020745", "NM_001605", "NM_005763"]
