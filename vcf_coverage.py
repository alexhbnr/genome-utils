#!/usr/bin/env python3

import argparse
import os.path
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
import vcf


def count_depths(vcf_file, *, include=[], exclude=[]):
    '''Walk over the VCF file and count the depths of all sites in all
    specified samples.
    '''
    reader = vcf.Reader(open(vcf_file, 'rb'))

    # determine which samples from a VCF to analyze
    include = set(reader.samples) if not include else set(include)
    exclude = set(exclude)
    sample_ids = include - exclude

    depth_counts = {s: defaultdict(int) for s in sample_ids}
    for rec in reader:
        for s in sample_ids:
            call = rec.genotype(s)
            if call.called:
                depth_counts[s][call['DP']] += 1

    # filter out samples without any information about coverage
    # and convert the final dict into a data frame
    return pd.DataFrame.from_dict({k: v for k, v in depth_counts.items()
                                        if None not in v.keys()})


def plot_single_sample(sample_id, depth_counts):
    '''Plot coverage distribution of SNPs in a given sample.'''
    fig = plt.figure()

    depth_counts[sample_id].plot(linewidth=3, figsize=(13, 9), fontsize=15)
    plt.title('SNP coverage in sample: ' + sample_id, fontsize=20)
    plt.xlabel('coverage', fontsize=16)
    plt.ylabel('count', fontsize=16)

    return fig



def plot_all_samples(depth_counts):
    '''Plot coverage distribution of SNPs in a given sample.'''
    depth_counts.plot(linewidth=3, figsize=(13, 9), fontsize=15)
    plt.title('SNP coverage in all samples', fontsize=20)
    plt.xlabel('coverage', fontsize=16)
    plt.ylabel('count', fontsize=16)
    plt.ylim(0, max(depth_counts.max()))
    print(max(depth_counts.max()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Plot coverage distribution of samples '
                                     'in a given VCF file')
    parser.add_argument('--vcf', help='VCF file to analyze', required=True)
    parser.add_argument('--separate', help='Plot distribution for each sample '
                        'separately', action='store_true')
    parser.add_argument('--include', help='Samples from VCF to analyze',
                        nargs='*', default=[])
    parser.add_argument('--exclude', help='Samples from VCF to skip',
                        nargs='*', default=[])
    args = parser.parse_args()

    depth_counts = count_depths(args.vcf,
                                include=args.include,
                                exclude=args.exclude)
    if args.separate:
        with PdfPages(os.path.basename(args.vcf) + '.pdf') as pdf:
            for sample in depth_counts.columns:
                plot_single_sample(sample, depth_counts)
                pdf.savefig()
                plt.close()
    else:
        plot_all_samples(depth_counts)
        plt.savefig(os.path.basename(args.vcf) + '.svg')

