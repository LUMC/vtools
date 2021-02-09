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

from vtools.gcoverage import fractions_at_least, qualmean

import numpy as np


def test_fractions_at_least():
    # The numbers are randomly generated but are are sorted for convenience:
    # list(sorted([random.randint(0, 100) for _ in range(100)]))
    hundred_numbers = [
        1, 2, 7, 8, 8, 9, 9, 9, 9, 11, 11, 14, 17, 19, 22, 23, 23, 24, 24, 25,
        26, 27, 27, 27, 27, 29, 29, 32, 34, 34, 35, 36, 37, 39, 39, 41, 41, 44,
        44, 44, 46, 47, 47, 48, 48, 50, 50, 51, 53, 54, 54, 56, 56, 56, 57, 58,
        59, 59, 63, 63, 65, 66, 67, 67, 67, 70, 70, 71, 73, 73, 76, 76, 77, 77,
        79, 79, 80, 81, 82, 84, 84, 85, 86, 86, 86, 91, 91, 93, 93, 93, 93, 94,
        94, 95, 95, 95, 96, 96, 97, 97]

    fractions = fractions_at_least(np.array(hundred_numbers, dtype=np.uint8),
                                   values_at_least=(10, 20, 30, 40, 50))
    assert fractions[0] == 91 / 100
    assert fractions[1] == 86 / 100
    assert fractions[2] == 73 / 100
    assert fractions[3] == 65 / 100
    assert fractions[4] == 55 / 100


def test_qualmean():
    qualities = np.array([5, 15, 25, 35, 45, 75, 95])
    assert math.isclose(qualmean(qualities), 12.9934485294235)


