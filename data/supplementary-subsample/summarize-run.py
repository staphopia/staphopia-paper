#! /usr/bin/env python3
# Summarize Staphopia run.
from collections import OrderedDict
import gzip
import json
SUMMARY = OrderedDict()
PROJECT = 'michellesu/read-lab'

def format_percent(number):
    return '{:.2f}'.format(number)


def milliseconds_to_minutes(total_time):
    return '{:.2f}'.format(total_time/1000/60)


def get_details():
    try:
        for job in api.tasks.query(project=PROJECT, status="COMPLETED").all():
            if job.name == SUMMARY['sample']:
                runtime = milliseconds_to_minutes(
                    job.execution_status.running_duration
                )
                SUMMARY['runtime'] = runtime
                price = str(job.price).split('amount=')
                if len(price) == 2:
                    price = price[1].replace('>', '')
                else:
                    price = "0.00"
                SUMMARY['price'] = price
    except SbgError:
        print('Unable to retrieve task info.')


def get_file_path(file, directory, sample):
    if '{1}' in file:
        file = file.format(directory, sample)
    else:
        file = file.format(directory)

    return os.path.abspath(os.path.realpath(file))


def parse_jellyfish(jf):
    """Get md5sum of a file."""
    from subprocess import Popen, PIPE
    f = Popen(['jellyfish', 'dump', '-c',  jf], stdout=PIPE)
    stdout, stderr = f.communicate()

    total = 0
    singleton = 0
    for line in stdout.decode("utf-8").split('\n'):
        if line:
            kmer, count = line.split()
            if int(count) == 1:
                singleton += 1
            total += 1

    SUMMARY['total_kmer'] = total
    SUMMARY['total_singleton'] = singleton


def parse_assembly(assembly):
    with open(assembly, 'r') as file_handle:
        json_data = json.load(file_handle)

    SUMMARY['total_contig'] = json_data['total_contig']
    SUMMARY['total_contig_length'] = json_data['total_contig_length']
    SUMMARY['max_contig_length'] = json_data['max_contig_length']
    SUMMARY['mean_contig_length'] = json_data['mean_contig_length']
    SUMMARY['median_contig_length'] = json_data['median_contig_length']
    SUMMARY['min_contig_length'] = json_data['min_contig_length']
    SUMMARY['n50_contig_length'] = json_data['n50_contig_length']


def parse_fastq(fastq, original=True):
    with open(fastq, 'r') as file_handle:
        json_data = json.load(file_handle)

    if original:
        SUMMARY['coverage_original'] = json_data['qc_stats']['coverage']
    else:
        SUMMARY['coverage_cleanup'] = json_data['qc_stats']['coverage']


def parse_annotation(annotation):
    with open(annotation, 'r') as file_handle:
        for line in file_handle:
            line = line.rstrip()
            if line.startswith('contigs'):
                SUMMARY['total_contig_200bp'] = line.replace('contigs: ', '')
            elif line.startswith('gene'):
                SUMMARY['total_gene'] = line.replace('gene: ', '')


def parse_variants(variant):
    total_snps = 0
    total_indel = 0
    with gzip.open(variant, 'r') as file_handle:
        for line in file_handle:
            line = line.decode("utf-8")
            if line.startswith("#"):
                continue
            cols = line.split('\t')
            ref = cols[3]
            alt = cols[4]

            if len(ref) == 1 and len(alt) == 1:
                total_snps += 1
            else:
                total_indel += 1

    SUMMARY['total_variant'] = total_snps + total_indel
    SUMMARY['total_snps'] = total_snps
    SUMMARY['total_indel'] = total_indel


if __name__ == '__main__':
    import os
    import argparse as ap
    import sevenbridges as sbg
    from sevenbridges.errors import SbgError

    parser = ap.ArgumentParser(
        prog='summarize-run.py',
        conflict_handler='resolve',
        description=('Summarize Staphopia run.'))

    parser.add_argument('directory', metavar="INPUT_DIRECTORY", type=str,
                        help='Directory to to read input files from.')
    args = parser.parse_args()

    sample = args.directory.split("cgc/")[1].replace('/', '')
    mutations = sample.split('-')[1]
    simulation = sample.split('-')[2]
    coverage = sample.split('-')[3]
    SUMMARY['sample'] = sample
    api = sbg.Api(config=sbg.Config(profile='cgc'))
    get_details()
    SUMMARY['mutations'] = mutations
    SUMMARY['simulation'] = simulation
    SUMMARY['coverage'] = coverage
    files = {
        'annotation_txt': '{0}/analyses/annotation/{1}.txt',
        'assembly_contig_stats': '{0}/analyses/assembly/{1}.contigs.json',
        'fastq_cleanup': '{0}/analyses/illumina-cleanup/{1}.cleanup.fastq.json',
        'fastq_original': '{0}/analyses/illumina-cleanup/{1}.original.fastq.json',
        'kmer': '{0}/analyses/kmer/{1}.jf',
        'variants': '{0}/analyses/variants/{1}.variants.vcf.gz',
    }

    order = ['annotation_txt', 'assembly_contig_stats', 'fastq_cleanup',
             'fastq_original', 'kmer', 'variants']

    for f in order:
        input_file = get_file_path(files[f], args.directory, sample)
        if f == 'kmer':
            parse_jellyfish(input_file)
        elif f == 'annotation_txt':
            parse_annotation(input_file)
        elif f == 'assembly_contig_stats':
            parse_assembly(input_file)
        elif f == 'fastq_cleanup':
            parse_fastq(input_file, original=False)
        elif f == 'fastq_original':
            parse_fastq(input_file)
        elif f == 'variants':
            parse_variants(input_file)

    header = []
    row = []
    for key, val in SUMMARY.items():
        header.append(key)
        row.append(str(val))

    print('\t'.join(header))
    print('\t'.join(row))
