#!/usr/bin/env python3

import argparse
import os.path
from collections import defaultdict

import vcf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

from bam_plotdamage import mismatch, MISMATCHES


def count_snps(vcf_file, sample, chrom):
    '''Analyze counts of SNPs in a sample in a given VCF file.'''
    reader = vcf.Reader(open(vcf_file, 'rb'))

    # initialize counters of mismatches from 5' and 3' ends of reads
    snp_counts = pd.DataFrame({
        'subst': MISMATCHES,
        'count': np.zeros(len(MISMATCHES), dtype=int)
    })

    for rec in reader.fetch(chrom):
        call = rec.genotype(sample)
        if not (rec.is_snp and call.called): continue

        for b in list(call.gt_bases):
            if rec.REF != b:
                key = mismatch(rec.REF, b)
                snp_counts.loc[snp_counts['subst'] == key, 'count'] += 1

    return snp_counts


def plot_snp_counts(snp_counts, sample, asfreq=False):
    '''Plot distribution of substitution frequencies.'''
    snp_counts['freq'] = snp_counts['count'] / sum(snp_counts['count'])

    which = 'freq' if asfreq else 'count'

    fig = plt.figure()
    sns.barplot(x='subst', y=which, data=snp_counts)
    plt.title('Distribution of SNPs in sample ' + sample)
    plt.xlabel('SNP type')
    plt.ylabel(which)

    return fig
    

def save_snp_counts(snp_counts, vcf_file, sample, asfreq):
    '''Save SNP count distribution plot to a file.'''
    fig = plot_snp_counts(snp_counts, sample, asfreq)
    output_file = 'snp_counts__' + sample + '__' + os.path.basename(vcf_file) + '.svg'
    plt.savefig(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Plot distribution of SNP counts in '
                                     'a sample from a given VCF file')
    parser.add_argument('--vcf', help='VCF file to analyze', required=True)
    parser.add_argument('--sample', help='Which sample to analyze',
                        required=True)
    parser.add_argument('--chrom', help='Chromosome', default='Y')
    parser.add_argument('--asfreq', help='Plot as frequencies or counts?',
                        action='store_true')
    parser.add_argument('--dir', help='Where to put the plot?', default='.')
    args = parser.parse_args()

    snp_counts = count_snps(args.vcf, args.sample, args.chrom)
    save_snp_counts(snp_counts, args.vcf, args.sample, args.asfreq)


