"""
setup.py
~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from os.path import abspath, dirname, join

from setuptools import setup, find_packages, Extension

readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()

setup(
    name="v-tools",
    version="1.1.0",
    description="Various tools operating over VCF files",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    author="Sander Bollen, Redmar van den Berg",
    author_email="KG_bioinf@lumc.nl",
    url="https://git.lumc.nl/klinische-genetica/capture-lumc/vtools",
    license="MIT",
    package_dir={"": "src"},
    packages=find_packages("src"),
    package_data={
        'vtools': ['vtools/*.pyx']
    },
    python_requires=">=3.6",
    zip_safe=False,
    include_package_data=True,
    setup_requires=["cython", "numpy"],
    install_requires=[
        "click",
        "cyvcf2",
        "numpy",
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
    ext_modules=[Extension("vtools.optimized", ["src/vtools/optimized.pyx"])]
)
