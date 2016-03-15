import argparse
import os.path
from collections import defaultdict
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


bases = 'ACGT'


def mismatch_key(b1, b2):
    '''Create a string X-Y to be used as a key to a data frame.'''
    return '-'.join((b1, b2))


def mismatch_combinations():
    '''Generate all combination of mismatches to be used as keys
    to a DataFrame.
    '''
    mismatch_comb = [mismatch_key(b1, b2) for b1 in bases
                                          for b2 in bases if b1 != b2]
    return mismatch_comb


def count_mismatches_in_read(mismatch_table, ref_bases, read_bases, len_limit):
        '''Calculate the number of mismatches along a read and update
        the mismatch table. Do not analyze sites within a read which are
        further from the beginning than len_limit.
        '''
        for pos, (b_ref, b_read) in enumerate(zip(ref_bases, read_bases)):
            # skip the rest of the read if already beyond the limit
            if pos >= len_limit: break

            # skip the base if it's not A, C, G or T
            if b_ref not in bases or b_read not in bases: continue

            # if there's a mismatch on this position, increment the counter
            if b_ref != b_read:
                mismatch_table.ix[pos, mismatch_key(b_ref, b_read)] += 1


def analyze_bam(bam_path, len_limit=20):
    '''Count the number of occurences of different substitutions in all
    reads as well as the total number of reads and the number of unmapped
    reads.
    '''   
    # initialize counters of reads and mismatch-counting data frames
    n_fwd_reads, n_rev_reads = 0, 0
    fwd_mismatches = pd.DataFrame(0, columns=mismatch_combinations(), index=range(len_limit))
    rev_mismatches = pd.DataFrame(0, columns=mismatch_combinations(), index=range(len_limit))
    
    ref_genome = '/mnt/solexa/Genomes/hg19_evan/whole_genome.fa'

    with pysam.AlignmentFile(bam_path, 'rb') as bamf:
        fastaf = pysam.FastaFile(ref_genome)
        
        for read in bamf:
            ref_bases = list(fastaf.fetch(reference=bamf.getrname(read.reference_id),
                                          start=read.reference_start,
                                          end=read.reference_end))
            read_bases = list(read.query_alignment_sequence)
            mismatch_table = fwd_mismatches

            # if the read is in a reversed orientation, reverse the fetched
            # sequences and change reference to the table to be updated
            if read.is_reverse:
                read_bases = reversed(read_bases)
                ref_bases = reversed(ref_bases)
                mismatch_table = rev_mismatches

            # update the mismatch table based on data from this read
            count_mismatches_in_read(mismatch_table, ref_bases, read_bases, len_limit)
            
            if read.is_reverse:
                n_rev_reads += 1
            else:
                n_fwd_reads += 1
    
    
    fwd_mismatches = fwd_mismatches.div(n_fwd_reads) * 100
    rev_mismatches = rev_mismatches.div(n_rev_reads) * 100
    
    print('Total number of forward reads: {}'.format(n_fwd_reads))
    print('Total number of reverse reads: {}'.format(n_rev_reads))
    
    return fwd_mismatches, rev_mismatches[::-1]


def plot_mismatches(mismatch_table, title, ylim=30):
    '''Plot the counts of mismatches in a specified orientation.'''
    fig = plt.figure()

    mismatch_table.plot(linewidth=3, figsize=(13, 9), fontsize=15)
    plt.title(title, fontsize=20)
    plt.xlabel('position from the beginning of the read [bp]', fontsize=16)
    plt.ylabel('fraction of mismatch type at a position [%]', fontsize=16)
    plt.xticks(range(len(mismatch_table)))
    plt.ylim(0, ylim)

    return fig

    
def save_mismatches(mismatch_table, strand, output_dir, bam_file, len_limit, ylim=30):
    title = strand + ' strand ' + '-- ' + os.path.basename(bam_file)
    fig = plot_mismatches(mismatch_table, title, ylim)
    output_file = strand + '_damage_' + os.path.basename(bam_file) + '.png'
    plt.savefig(os.path.join(output_dir, output_file))


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Analyze damage patterns in a BAM file')
    parser.add_argument('--bam', help='BAM file to analyze', required=True)
    parser.add_argument('--len_limit', help='How deep into a read to look for damage?', default=20)
    parser.add_argument('--dir', help='Where to put the plot?', default='.')
    parser.add_argument('--format', help='Format of the output image',
                        default='png')
    args = parser.parse_args()

    fwd_mismatches, rev_mismatches = analyze_bam(args.bam, args.len_limit)
    save_mismatches(fwd_mismatches, 'forward', args.dir, args.bam, args.len_limit)
    save_mismatches(rev_mismatches, 'reverse', args.dir, args.bam, args.len_limit)
