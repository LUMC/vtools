# Copyright (c) 2021 Leiden University Medical Center
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

import tempfile
import zipfile
from pathlib import Path

from vtools.util import zip_to_files


def test_zip_to_files():
    temp_dir = tempfile.mkdtemp()
    test_zip = Path(temp_dir, "test.zip")
    temp_file = Path(temp_dir, "test.txt")
    temp_file.write_text("Test\n")
    with zipfile.ZipFile(test_zip, "w") as testzip_h:
        for i in range(10):
            testzip_h.write(temp_file, "test" + str(i))
    count = 0
    for file in zip_to_files(test_zip):
        assert Path(file).read_text() == "Test\n"
        count += 1
    assert count == 10
