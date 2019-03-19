"""
setup.py
~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from os.path import abspath, dirname, join

from setuptools import setup, find_packages


try:
    from Cython.Build import cythonize
except ImportError:
    raise NotImplementedError("Installing cython on the fly not yet supported")


try:
    import numpy as np
except ImportError:
    raise NotImplementedError("Installing numpy on the fly not yet supported")

readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()

# create extensions and add numpy includes to all of them.
cython_extensions = cythonize("vtools/*.pyx")
for ext in cython_extensions:
    ext.include_dirs.append(np.get_include())

setup(
    name="vtools",
    version="0.0.1",
    description="Various tools operating over VCF files",
    author="Sander Bollen",
    author_email="a.h.b.bollen@lumc.nl",
    license="MIT",
    packages=find_packages(),
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
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    ext_modules=cython_extensions
)
