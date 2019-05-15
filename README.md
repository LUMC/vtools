[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/vtools/README.html)

vtools
======

Little toolset operating over VCF files. Uses cyvcf2 and cython under
the hood for speed.

## Installation

### PyPI
vtools is now on pip! Since the 'vtools' name is already taken by another
package, installing _this_ vtools requires installing the following:

```bash
pip install v-tools
```

After installation, tools will still be called `vtools-<tool>`. Programmatic
access also simply works with

```python
import vtools
```

### Conda

```bash
conda install -c bioconda vtools
```


Tools
-----

### vtools-filter

Filter VCF files based on a few criteria. Will output both a filtered VCF
file, and a VCF file containing all the filtered-out variants.

####  Filter criteria

| name | meaning | optional |
| ---- | ------- | -------- |
| NON_CANONICAL | Non-canonical chromosome | Yes |
| INDEX_UNCALLED | Index uncalled or homozygous reference | Yes |
| TOO_HIGH_GONL_AF | Too high GonL allele frequency | Yes |
| TOO_HIGH_GNOMAD_AF | Too high GnomAD allele frequency | Yes |
| LOW_GQ | Too low GQ on index sample | Yes |
| DELETED_ALLELE | The only ALT allele is a deleted allele | No |

#### Configuration 

Configuration of filters goes by a little JSON file. 
See [here](cfg/example-filter.json) for an example.


#### Usage

```bash
Usage: vtools-filter [OPTIONS]

Options:
  -i, --input PATH                Path to input VCF file  [required]
  -o, --output PATH               Path to output (filtered) VCF file
                                  [required]
  -t, --trash PATH                Path to trash VCF file  [required]
  -p, --params-file PATH          Path to filter params json  [required]
  --index-sample TEXT             Name of index sample  [required]
  --immediate-return / --no-immediate-return
                                  Immediately write filters to file upon
                                  hitting one filter criterium. Default = True
  --help                          Show this message and exit.

```

### vtools-stats

Collects some general statistics about a VCF file, and writes a json to
stdout.

#### Usage

```bash
Usage: vtools-stats [OPTIONS]

Options:
  -i, --input FILE  Input VCF file  [required]
  --help            Show this message and exit.
```

### vtools-gcoverage

Collect coverage metrics over a gVCF file for every exon or every transcript
in a refFlat file. This assumes the input VCF file is at least similar to
GATK's gVCF files. gVCF files are only expected to have one sample; if
your input file contains multiple samples, we simply take the first only.

Output is a simple TSV file with the following columns

| column | meaning |
| ------ | ------- |
| exon | exon number |
| gene | gene name / symbol / id |
| mean_dp | mean DP value over the exon |
| mean_gq | mean GQ value over the exon* |
| median_dp | median DP value over the exon |
| median_gq | median GQ value over the exon |
| perc_at_least_{10, 20, 30, 50, 100}_dp | Percentage of exon with DP value over value |
| perc_at_least_{10, 29, 30, 50, 90}_gq | Percentage of exon with GQ value over exon | 
| transcript | transcript name / symbol / id |

*: mean GQ value is computed by first calculating the P-value of all GQ 
values, then calculating the mean over these P-values, and lastly 
converting this number back to a phred score.

#### Usage

```bash
Usage: vtools-gcoverage [OPTIONS]

Options:
  -I, --input-gvcf PATH          Path to input VCF file  [required]
  -R, --refflat-file PATH        Path to refFlat file  [required]
  --per-exon / --per-transcript  Collect metrics per exon or per transcript
  --help                         Show this message and exit.
```

### vtools-evaluate

Evaluate a VCF file to a baseline VCF file containing true positives. 
We only consider variants that are present in both VCF files. This makes
it useful when the two VCF files have been produced by wildly different
technologies. E.g, when comparing a WES VCF file vs a SNP array, this
tool can be quite useful.

Output is a simple JSON file listing counts of concordant and discordant
alleles. 

Multisample VCF files are allowed; the samples to be evaluated have to be set 
through a CLI argument.


#### Usage

```bash
Usage: vtools-evaluate [OPTIONS]

Options:
  -c, --call-vcf PATH           Path to VCF with calls to be evaluated
                                [required]
  -p, --positive-vcf PATH       Path to VCF with known calls  [required]
  -cs, --call-samples TEXT      Sample(s) in call-vcf to consider. May be
                                called multiple times  [required]
  -ps, --positive-samples TEXT  Sample(s) in positive-vcf to consider. May be
                                called multiple times  [required]
  --help                        Show this message and exit.
```

## License

MIT
