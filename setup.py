"""
setup.py
~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from os.path import abspath, dirname, join
import sys
import pkg_resources
import subprocess

from setuptools import setup, find_packages


# Temporarily install dependencies required by setup.py before trying to
# import them. From https://bitbucket.org/dholth/setup-requires
sys.path[0:0] = ['setup-requires']
pkg_resources.working_set.add_entry('setup-requires')


def missing_requirements(specifiers):
    for specifier in specifiers:
        try:
            pkg_resources.require(specifier)
        except pkg_resources.DistributionNotFound:
            yield specifier


def install_requirements(specifiers):
    to_install = list(specifiers)
    if to_install:
        cmd = [sys.executable, "-m", "pip", "install",
               "-t", "setup-requires"] + to_install
        subprocess.call(cmd)


requires = ['cython', 'numpy']
install_requirements(missing_requirements(requires))


from Cython.Build import cythonize
import numpy as np

readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()

# create extensions and add numpy includes to all of them.
cython_extensions = cythonize("vtools/*.pyx")
for ext in cython_extensions:
    ext.include_dirs.append(np.get_include())

setup(
    name="v-tools",
    version="1.0.0",
    description="Various tools operating over VCF files",
    long_description=long_desc,
    author="Sander Bollen",
    author_email="a.h.b.bollen@lumc.nl",
    url="https://git.lumc.nl/klinische-genetica/capture-lumc/vtools",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.6",
    zip_safe=False,
    include_package_data=True,
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
            "vtools-gcoverage = vtools.cli:gcoverage_cli",
            "vtools-evaluate = vtools.cli:evaluate_cli"
        ]
    },
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    ext_modules=cython_extensions
)
