vtools
======

Little toolset operating over VCF files. Uses cyvcf2 and cython under
the hood for speed.


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

Configuration of filters goes by a little JSON file. See [here]() for an 
example, and [here]() for the json schema.


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
  
