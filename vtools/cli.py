"""
vtools.cli
~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) Leiden University Medical Center
:license: MIT
"""
import json
import click
from cyvcf2 import VCF, Writer

from .evaluate import site_concordancy
from .filter import FilterParams, FilterClass, Filterer
from .stats import Stats
from .gcoverage import RefRecord, region_coverages


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
def evaluate_cli(call_vcf, positive_vcf, call_samples, positive_samples):
    c_vcf = VCF(call_vcf, gts012=True)
    p_vcf = VCF(positive_vcf, gts012=True)
    evaluated = site_concordancy(c_vcf, p_vcf, call_samples, positive_samples)
    print(json.dumps(evaluated))


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
@click.option("-I", "--input-gvcf",
              type=click.Path(exists=True, readable=True),
              required=True,
              help="Path to input VCF file")
@click.option("-R", "--refflat-file",
              type=click.Path(exists=True, readable=True),
              required=True,
              help="Path to refFlat file")
@click.option("--per-exon/--per-transcript",
              default=True,
              help="Collect metrics per exon or per transcript")
def gcoverage_cli(input_gvcf, refflat_file, per_exon):
    reader = VCF(input_gvcf)
    header = None
    with open(refflat_file) as handle:
        for line in handle:
            r = RefRecord.from_line(line)
            if not per_exon:
                regions = [x[1] for x in r.cds_exons]
                cov = region_coverages(reader, regions)
                cov['transcript'] = r.transcript
                cov['gene'] = r.gene
                if header is None:
                    header = "\t".join(sorted(cov.keys()))
                    print(header)
                vals = [str(cov[k]) for k in sorted(cov.keys())]
                print("\t".join(vals))
            else:
                for i, reg in enumerate(r.exons):
                    cov = region_coverages(reader, [reg])
                    cov['transcript'] = r.transcript
                    cov['gene'] = r.gene
                    if r.forward:
                        cov['exon'] = i+1
                    else:
                        cov['exon'] = len(r.exons) - i
                    if header is None:
                        header = "\t".join(sorted(cov.keys()))
                        print(header)
                    vals = [str(cov[k]) for k in sorted(cov.keys())]
                    print("\t".join(vals))
