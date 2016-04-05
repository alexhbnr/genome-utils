#!/usr/bin/env python3

import sys
import argparse
import random
from itertools import islice, chain
from os.path import basename, splitext

import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def bquals_in_read(read):
    '''Return a list of base qualities in a given read.'''
    return [ord(q) - 33 for q in read.qual]


def subsample_reads(bam, sample_size):
    '''Perform reservoir sampling to get a subset of reads from a BAM file.'''
    # take the first sample of candidate reads of size 'sample_size'
    sampled_reads = [read for read in islice(bam.fetch(), sample_size)]

    for i, read in enumerate(bam.fetch()):
        # skip reads that are already part of the initial sample
        if i < sample_size: continue
        
        # include each following i-th read with a probability of k / i
        k = random.randint(0, i - 1)
        if k < sample_size:
            sampled_reads[k] = read

    return sampled_reads


def bquals_sample(bamfile, n_reads):
    '''Get a sample of base qualities from a given BAM file.'''
    bam = pysam.AlignmentFile(bamfile)
    reads = subsample_reads(bam, n_reads)
    bquals = np.fromiter(chain.from_iterable(bquals_in_read(r) for r in reads),
                         np.int)
    return bquals


def plot_bqual_dist(bquals, bamfile):
    '''Plot the distribution of base qualities.'''
    max_bqual = 60

    fig = plt.figure(figsize=(13, 9))
    plt.hist(bquals, bins=max_bqual, range=(0, max_bqual), normed=True)
    plt.title('Distribution of base qualities in ' + basename(bamfile), fontsize=20)
    plt.xlabel('base quality [Phred scaled]', fontsize=16)
    plt.ylabel('probability density', fontsize=16)

    return fig


def save_bqual_dist(bquals, bamfile):
    '''Save distribution of base qualities to a file.'''
    fig = plot_bqual_dist(bquals, bamfile)
    plt.savefig('bqdist_' + basename(bamfile) + '.svg', format='svg')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze the distribution of '
        'base qualities in a subset of reads from a given BAM file')
    parser.add_argument('--bam', help='BAM file to analyze', required=True)
    parser.add_argument('--subsample', help='Number of reads to sample',
                        default=10000)
    args = parser.parse_args()

    bquals = bquals_sample(args.bam, args.subsample)
    save_bqual_dist(bquals, args.bam)
