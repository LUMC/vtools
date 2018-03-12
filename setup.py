"""
setup.py
~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

from os.path import abspath, dirname, join

from setuptools import setup
from Cython.Build import cythonize

from vtools import __version__


readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()


setup(
    name="vtools",
    version=__version__,
    description="Various tools operating over VCF files",
    author="Sander Bollen",
    author_email="a.h.b.bollen@lumc.nl",
    license="MIT",
    packages=["vtools"],
    install_requires=[
        "click",
        "cyvcf2",
        "numpy",
        "cython",
        "tqdm"
    ],
    entry_points={
        "console_scripts": [
            "vtools-filter = vtools.cli:filter_cli",
            "vtools-stats = vtools.cli:stats_cli",
            "vtools-gcoverage = vtools.cli:gcoverage_cli"
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    ext_modules=cythonize("vtools/*.pyx")
)