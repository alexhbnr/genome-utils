#!/usr/bin/env python3

import argparse
import os.path
from collections import defaultdict

import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def count_read_lengths(bam_path):
    '''Count the length distribution of reads in a specified BAM file.'''
    length_counter = defaultdict(int)

    with pysam.AlignmentFile(bam_path, 'rb') as bamf:
        for read in bamf:
            length_counter[read.query_length] += 1

    return pd.Series(length_counter)


def plot_length_dist(read_lengths, title):
    '''Plot distribtion of read lengths and return the figure object
    (this allows to use this function in an interactive session).
    '''
    fig = plt.figure()

    read_lengths.plot(linewidth=3, figsize=(10, 6), fontsize=15)
    plt.title(title, fontsize=20)
    plt.xlabel('read length [bp]', fontsize=16)
    plt.ylabel('count', fontsize=16)
    plt.xlim(0, max(read_lengths.index))

    return fig


def save_length_dist(read_lengths, output_dir, bam_file):
    '''Plot the read length distribution to a file.'''
    fig = plot_length_dist(read_lengths, os.path.basename(bam_file))

    output_img = 'read_lengths_' + os.path.basename(bam_file) + '.svg'
    plt.savefig(os.path.join(output_dir, output_img))


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Calculate distribution of read lengths '
                                     'in the given BAM file.')
    parser.add_argument('--bam', help='BAM file to analyze', required=True)
    parser.add_argument('--dir', help='Where to put the plot?', default='.')
    args = parser.parse_args()

    read_lengths = count_read_lengths(args.bam)
    save_length_dist(read_lengths, args.dir, args.bam)
