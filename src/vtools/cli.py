"""
vtools.cli
~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) Leiden University Medical Center
:license: MIT
"""
import gzip
import json
import os

import click

from cyvcf2 import VCF, Writer  # type: ignore

from .evaluate import site_concordancy
from .filter import FilterClass, FilterParams,  Filterer
from .gcoverage import refflat_and_gvcfs_to_tsv
from .stats import Stats


@click.command()
@click.option("-c", "--call-vcf", type=click.Path(exists=True),
              help="Path to VCF with calls to be evaluated",
              required=True)
@click.option("-p", "--positive-vcf", type=click.Path(exists=True),
              help="Path to VCF with known calls",
              required=True)
@click.option("-cs", "--call-samples", type=click.STRING, multiple=True,
              help="Sample(s) in call-vcf to consider. "
                   "May be called multiple times",
              required=True)
@click.option("-ps", "--positive-samples", type=click.STRING, multiple=True,
              help="Sample(s) in positive-vcf to consider. "
                   "May be called multiple times",
              required=True)
@click.option("-s", "--stats", type=click.Path(writable=True),
              help="Path to output stats json file")
@click.option("-dvcf", "--discordant-vcf", type=click.Path(writable=True),
              help="Path to output the discordant vcf file",
              required=False)
@click.option("-mq", "--min-qual", type=float,
              help="Minimum quality of variants to consider", default=30)
@click.option("-md", "--min-depth", type=int,
              help="Minimum depth of variants to consider", default=0)
def evaluate_cli(call_vcf, positive_vcf, call_samples, positive_samples,
                 min_qual, min_depth, stats, discordant_vcf):
    c_vcf = VCF(call_vcf, gts012=True)
    p_vcf = VCF(positive_vcf, gts012=True)
    st, disc = site_concordancy(c_vcf, p_vcf, call_samples,
                                positive_samples, min_qual, min_depth)
    # Write the stats json file
    if stats is None:
        print(json.dumps(st))
    else:
        with click.open_file(stats, 'w') as fout:
            fout.write(json.dumps(st))

    # If there were discordand records, and a discordant-vcf should be written
    if len(disc) > 0 and discordant_vcf:
        # make sure the parent folder exists
        parent_folder = os.path.dirname(discordant_vcf)
        os.makedirs(parent_folder, exist_ok=True)

        with click.open_file(discordant_vcf, 'w') as fout:
            # First, we write the vcf header
            with gzip.open(call_vcf, 'rt') as fin:
                for line in fin:
                    if line.startswith('#'):
                        fout.write(line)
                    else:
                        break
            # Then we write the vcf records that were discordant
            for record in disc:
                fout.write(str(record))


@click.command()
@click.option("-i", "--input", type=click.Path(exists=True),
              help="Path to input VCF file", required=True)
@click.option("-o", "--output", type=click.Path(writable=True),
              help="Path to output (filtered) VCF file", required=True)
@click.option("-t", "--trash", type=click.Path(writable=True),
              help="Path to trash VCF file", required=True)
@click.option("-p", "--params-file", type=click.Path(exists=True),
              help="Path to filter params json", required=True)
@click.option('--index-sample', type=click.STRING,
              help="Name of index sample", required=True)
@click.option("--immediate-return/--no-immediate-return",
              default=True,
              help="Immediately write filters to file "
                   "upon hitting one filter criterium. "
                   "Default = True")
def filter_cli(input, output, trash, params_file,
               index_sample, immediate_return):
    vcf = VCF(input, gts012=True)

    idx = vcf.samples.index(index_sample)
    for filter_item in list(FilterClass):
        vcf.add_filter_to_header(filter_item.value)

    out = Writer(output, vcf)
    tr = Writer(trash, vcf)

    filter_params = FilterParams(params_file)

    filter_it = Filterer(vcf, filter_params, idx, immediate_return)

    for record, fi in filter_it:
        if fi is None or len(fi) == 0:
            out.write_record(record)
        else:
            record.FILTER = [x.name for x in fi]
            tr.write_record(record)

    out.close()
    tr.close()


@click.command()
@click.option("-i",
              "--input",
              type=click.Path(exists=True, dir_okay=False, readable=True),
              required=True,
              help="Input VCF file")
def stats_cli(input):
    stats = Stats(input)
    print(stats.as_json)


@click.command()
@click.argument("input-gvcf",
                nargs=-1,
                type=click.Path(exists=True, readable=True),
                required=True)
@click.option("-R", "--refflat-file",
              type=click.Path(exists=True, readable=True),
              required=True,
              help="Path to refFlat file")
@click.option("--per-exon/--per-transcript",
              default=True,
              help="Collect metrics per exon or per transcript")
def gcoverage_cli(input_gvcf, refflat_file, per_exon):
    for line in refflat_and_gvcfs_to_tsv(refflat_file, input_gvcf, per_exon):
        print(line)
