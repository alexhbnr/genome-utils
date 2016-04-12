#!/usr/bin/env python3

import sys
import argparse
import os.path
import random
from itertools import islice
from collections import defaultdict

import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')


BASES = 'ACGT'

# translation dict for converting DNA sequence into reverse complement
TR_TABLE = str.maketrans(BASES, BASES[::-1])


def mismatch(b1, b2):
    '''Create a string "X-Y" to be used as a key to a data frame.'''
    return b1 + '-' + b2


# all possible combinations of mismatches
MISMATCHES = [mismatch(b1, b2) for b1 in BASES
                               for b2 in BASES if b1 != b2]


def mismatches_in_read(mism5_count, mism3_count, ref5_count, ref3_count,
                       ref_bases, read_bases, len_limit):
    '''Count the number of mismatches along a read as well as the
    composition of reference bases.
    '''
    read_len = len(read_bases)

    for pos, (b_ref, b_read) in enumerate(zip(ref_bases, read_bases)):
        rev_pos = read_len - pos - 1

        # skip the position if the reference does not have a valid base
        if b_ref not in BASES or b_read not in BASES: continue

        # increment counters of ref alleles at a position
        if pos     < len_limit: ref5_count.ix[pos,     b_ref] += 1
        if rev_pos < len_limit: ref3_count.ix[rev_pos, b_ref] += 1

        # if there's a mismatch on this position, increment the counter
        if b_ref != b_read:
            key = mismatch(b_ref, b_read)
            if pos     < len_limit: mism5_count.ix[pos, key] += 1
            if rev_pos < len_limit: mism3_count.ix[rev_pos, key] += 1


def calc_frequencies(mism_counts, ref_count):
    ''' Calculate substitution frequencies from counts of observed
    substitution patterns and counts of reference bases at each
    position.
    '''
    return pd.concat((mism_counts.filter(regex=r'^' + b).divide(ref_count[b], axis=0)
                      for b in BASES), axis=1)


def count_mismatches(bam_path, len_limit=30):
    '''Count the number of occurences of different substitutions in all
    reads as well as the total number of reads and the number of unmapped
    reads.
    '''   
    # initialize counters of mismatches from 5' and 3' ends of reads
    mism5_counts = pd.DataFrame(0, columns=MISMATCHES, index=range(len_limit))
    mism3_counts = pd.DataFrame(0, columns=MISMATCHES, index=range(len_limit))

    # initialize counters of reference bases from the 5' and 3' ends of reads
    ref5_counts = pd.DataFrame(0, columns=list(BASES), index=range(len_limit))
    ref3_counts = pd.DataFrame(0, columns=list(BASES), index=range(len_limit))
    
    ref_genome = '/mnt/solexa/Genomes/hg19_evan/whole_genome.fa'

    with pysam.AlignmentFile(bam_path, 'rb') as bamf:
        fastaf = pysam.FastaFile(ref_genome)

        sample_size = 100000
        sampled_reads = [read for read in islice(bamf, sample_size)]
        
        for (i, read) in enumerate(bamf):
            if i < sample_size: continue

            k = random.randint(0, i - 1)
            if k < sample_size:
                sampled_reads[k] = read

        print('Subsampling of ' + str(sample_size) + ' reads finished. '
              'Analyzing damage patterns...')

        for (i, read) in enumerate(sampled_reads):
            if i % 5000 == 0: print(i, ' reads analyzed...', end='\r')
            ref_bases = fastaf.fetch(read.reference_name,
                                     read.reference_start,
                                     read.reference_end)
            read_bases = read.query_alignment_sequence

            if read.is_reverse:
                ref_bases = ref_bases[::-1].translate(TR_TABLE)
                read_bases = read_bases[::-1].translate(TR_TABLE)

            mismatches_in_read(mism5_counts, mism3_counts,
                               ref5_counts, ref3_counts,
                               ref_bases, read_bases, len_limit)

    mism5_freqs = calc_frequencies(mism5_counts, ref5_counts)
    mism3_freqs = calc_frequencies(mism3_counts, ref3_counts)

    return mism5_freqs, mism3_freqs[::-1]


def plot_mismatches(mism_freqs, read_end, bam_file):
    '''Plot the frequencies of mismatches.'''
    fig = plt.figure()

    end_str = str(read_end) + '\' end'

    mism_freqs.plot(linewidth=3, figsize=(13, 9), fontsize=15)
    plt.title('mismatches from the ' + end_str + ' -- ' +
              os.path.basename(bam_file), fontsize=20)
    plt.xlabel('position from the ' + end_str, fontsize=16)
    plt.ylabel('proportion of mismatches', fontsize=16)
    plt.xticks(range(0, len(mism_freqs), 5))

    return fig

    
def save_mismatches(mismatch_table, read_end, output_dir, bam_file, len_limit):
    fig = plot_mismatches(mismatch_table, read_end, bam_file)

    end_str = str(read_end) + '_prime_end'

    output_file = 'subst_from_' + end_str + '__' + os.path.basename(bam_file) + '.svg'
    plt.savefig(os.path.join(output_dir, output_file))


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Analyze damage patterns in a BAM file')
    parser.add_argument('--bam', help='BAM file to analyze', required=True)
    parser.add_argument('--len_limit', help='How deep into reads to look '
                        'for damage?', type=int, default=30)
    parser.add_argument('--which', help='Which substitutions to plot?',
                        nargs='*', default=MISMATCHES)
    parser.add_argument('--dir', help='Where to put the plot?', default='.')
    args = parser.parse_args()

    if not set(args.which).issubset(MISMATCHES):
        print('The only valid substitution patterns are:', ', '.join(MISMATCHES))
        sys.exit()

    mism5_freqs, mism3_freqs = count_mismatches(args.bam, args.len_limit)
    save_mismatches(mism5_freqs[args.which], 5, args.dir, args.bam, args.len_limit)
    save_mismatches(mism3_freqs[args.which], 3, args.dir, args.bam, args.len_limit)
