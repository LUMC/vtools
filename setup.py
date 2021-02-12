"""
setup.py
~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from pathlib import Path

import pkg_resources

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

long_desc = (Path(__file__).parent / "README.md").read_text()


class BuildExt(build_ext):
    def build_extension(self, ext: Extension):
        # Only get numpy location inside setup_requires environment during
        # the build_ext phase.
        numpy_include = str(
            Path(pkg_resources.get_distribution("numpy").location,
                 "numpy", "core", "include"))
        ext.include_dirs = [numpy_include]
        super().build_extension(ext)


setup(
    name="v-tools",
    version="1.1.1-dev",
    description="Various tools operating over VCF files",
    long_description=long_desc,
    cmdclass={"build_ext": BuildExt},
    long_description_content_type="text/markdown",
    author="Sander Bollen, Redmar van den Berg",
    author_email="KG_bioinf@lumc.nl",
    url="https://github.com/lumc/vtools",
    license="MIT",
    package_dir={"": "src"},
    packages=find_packages("src"),
    package_data={
        'vtools': ['vtools/*.pyx']
    },
    python_requires=">=3.6",
    zip_safe=False,
    include_package_data=True,
    install_requires=[
        "click",
        "cyvcf2",
        "numpy>=1.19",
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
