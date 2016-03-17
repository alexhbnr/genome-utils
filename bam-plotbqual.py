import sys
from os.path import basename, splitext
import argparse
import random
import pysam
from itertools import islice, chain
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def get_bquals(read):
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


def plot_dist(bquals, bamfile, plotfile):
    '''Plot the distribution of base qualities.'''
    max_bqual = 60

    plt.hist(bquals, bins=max_bqual, range=(0, max_bqual), normed=True)
    plt.title('Distribution of base qualities in ' + basename(bamfile))
    plt.xlabel('base quality [Phred scaled]')
    plt.ylabel('probability density')

    # determine the plot output format from the specified image extension
    plot_format = splitext(plotfile)[1][1:]

    plt.savefig(plotfile, format=plot_format)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Analyze the distribution of '
        'base qualities in a subset of reads from a given BAM file')
    parser.add_argument('--bam', help='BAM file to analyze', required=True)
    parser.add_argument('--plotfile', help='Name of the output plot file '
                        '(the format can be png, pdf, ps, eps, svg)',
                        required=True)
    parser.add_argument('--subsample', help='Number of reads to sample',
                        default=10000)
    args = parser.parse_args(argv if argv else sys.argv[1:])

    bam = pysam.AlignmentFile(args.bam)
    reads = subsample_reads(bam, args.subsample)

    bquals = np.fromiter(chain.from_iterable(get_bquals(r) for r in reads),
                         np.int)
    plot_dist(bquals, args.bam, args.plotfile)


if __name__ == '__main__':
    main()
