#!/usr/bin/env python3

import argparse
import random
import sys
import signal
import functools

import vcf
from pybedtools import BedTool

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def sample_bases(reader, sample_name, print_fn, chrom, start=None, end=None):
    '''Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    '''
    for record in reader.fetch(chrom, start, end):
        if record.is_snp or record.is_monomorphic:
            sample_call = record.genotype(sample_name)
            called_gt = int(random.choice(sample_call.gt_alleles))

            ref_base = record.REF
            alt_base = record.alleles[called_gt]

            print_fn(record.CHROM, record.POS, ref_base, alt_base, called_gt)


def sample_in_regions(reader, bed, sample_name, print_fn):
    '''Sample alleles from the BAM file at each position specified in a BED
    file. Return the result as a list of tuples in the form of
    (chromosome, position, ref_base, called_base).
    '''
    for region in bed:
        sample_bases(reader, sample_name, print_fn, region.chrom, region.start, region.end)


def print_vcf_header(sample_name, handle):
    '''Print VCF header.'''
    print('##fileformat=VCFv4.1\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}'.
          format(sample=sample_name), file=handle)


def print_record(chrom, pos, ref_base, called_base, gt, rec_fmt, handle):
    '''Print information about sampled site in a given string format.'''
    alt = '.' if ref_base == called_base else called_base
    print(rec_fmt.format(chrom=chrom, start=pos - 1, end=pos, pos=pos,
                         ref=ref_base, allele=called_base, alt=alt, gt=gt))


def main(argv=None):
    parser = argparse.ArgumentParser(description='Sample alleles from a given'
        ' VCF file.')
    parser.add_argument('--vcf', help='VCF file to sample from', required=True)
    parser.add_argument('--chrom', help='Chromosome to sample', required=True)
    parser.add_argument('--bed', help='BED file with coordinates of regions'
                        '/sites to sample')
    parser.add_argument('--output', help='Name of the output file '
                        '(direct output to stdout if missing)', default=None)
    parser.add_argument('--format', help='Output as VCF or BED?',
                        choices=['VCF', 'BED'], required=True)
    parser.add_argument('--sample-name', help='Sample to take from the VCF')

    # if there were no arguments supplied to the main function, use sys.argv
    # (skipping the first element, i.e. the name of this script)
    args = parser.parse_args(argv if argv else sys.argv[1:])

    if args.format == 'VCF' and not args.sample_name:
        parser.error('Sample has to be specified')

    reader = vcf.Reader(open(args.vcf, 'rb'))

    # output the results as specified by user
    handle = open(args.output, 'w') if args.output else sys.stdout

    if args.format == 'VCF':
        print_vcf_header(args.sample_name, handle)
        rec_fmt = '{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}'
    else:
        rec_fmt = '{chrom}\t{start}\t{end}\t{ref}\t{allele}'

    print_fn = functools.partial(print_record, rec_fmt=rec_fmt, handle=handle)

    # if user specified a BED file, perform pileup on each region in that file
    if args.bed:
        bed = BedTool(args.bed)
        sample_in_regions(reader, bed, args.sample_name, print_fn)
    else: # otherwise scan the whole BAM file directly
        sample_bases(reader, args.sample_name, print_fn, args.chrom)

    if args.output:
        handle.close()


if __name__ == "__main__":
    main()
