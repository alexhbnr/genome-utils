#!/usr/bin/env python3

import random
import argparse
import sys
from collections import Counter
import pysam
from pybedtools import BedTool


def terminal_base(pos_in_read, read_len, reverse_strand, library_prep):
    '''Check if a position in a read is on a terminal end of the read
    that is relevant for filtering based on a specified library preparation
    method.'''
    is_terminal = False
    if library_prep == 'USER':
        if not reverse_strand and ((pos_in_read == 0) or (read_len - pos_in_read <= 2)): is_terminal = True
        if     reverse_strand and ((pos_in_read  < 2) or (read_len - pos_in_read == 1)): is_terminal = True
    if library_prep == 'non-USER_term3':
        if (pos_in_read < 3) or (read_len - pos_in_read <= 3): is_terminal = True
    return is_terminal


def damage_at_site(pileup_site, library_prep):
    '''Ignore the base in the pileup of reads in case that:
       A) - reference is C
          - and: 
              - USER treated: read shows T at the first 5' or last two 3' positions
              - non-USER treated:
                    a) read shows T at first three or last three positions
                    b) read shows T anywhere on the forward read
       B) - reference is G
          - and:
              - USER treated: read shows A at the first two 5' or last 3' positions
              - non-USER treated:
                    a) read shows A at first three or last three positions
                    b) read shows A anywhere on the reverse read
    '''
    ref_base, read_base, pos_in_read, read_len, reverse_strand = pileup_site

    damaged = False
    if library_prep in ['USER', 'non-USER_term3'] and  (ref_base == 'C' and read_base == 'T' and not reverse_strand and terminal_base(pos_in_read, read_len, reverse_strand, library_prep)): damaged = True
    if library_prep in ['USER', 'non-USER_term3'] and  (ref_base == 'G' and read_base == 'A' and     reverse_strand and terminal_base(pos_in_read, read_len, reverse_strand, library_prep)): damaged = True
    if library_prep == 'non-USER_all'             and ((ref_base == 'C' and read_base == 'T' and not reverse_strand) \
                                                    or (ref_base == 'G' and read_base == 'A' and     reverse_strand)): damage_state = True
    return damaged

 
def filter_out_damage(pileup_column, library_prep):
    '''Filter out bases in a given pileup column that are likely result
    of DNA damage.
    '''
    return [site for site in pileup_column if not damage_at_site(site, library_prep)]


def get_pileup_info(chrom, start, end, bam, ref_base):
    '''Return a reference base and a list of bases at a given site.'''
    site_pileup = []

    for col in bam.pileup(chrom, start, end):
        # analyze only columns within a specified region (pysam performs
        # pileup on whole read lengths overlapping a given region otherwise)
        if not (start <= col.pos < end): continue

        # walk through all reads overlapping the current column
        # and accumulate bases at that position
        for pileup_read in col.pileups:
            # skip deletions
            if pileup_read.is_del: continue

            pos_in_read = pileup_read.query_position
            read_len = pileup_read.alignment.query_length
            read_base = pileup_read.alignment.query_sequence[pos_in_read]
            is_reverse = pileup_read.alignment.is_reverse

            if read_base in "ACGT": 
                site_pileup.append((ref_base,
                                    read_base,
                                    pos_in_read,
                                    read_len,
                                    is_reverse))

    return site_pileup


def call_base(pileup_info, method):
    """Return the most frequently occuring element of a list."""
    bases = [base for _, base, _, _, _ in pileup_info]

    if method == 'majority':
        counts = Counter(bases).most_common()
        # take all bases with the highest count
        max_freq = max(c[1] for c in counts)
        bases = [c[0] for c in counts if c[1] == max_freq]

    return random.choice(bases)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Sample bases from BAM file')
    parser.add_argument('--bam', help='BAM file to sample from', required=True)
    parser.add_argument('--bed', help='BED file with sites to sample at',
                        required=True)
    parser.add_argument('--ref', help='FASTA reference', required=True)
    parser.add_argument('--output', help='Name of the output file '
                        '(write to stdout if missing)', default=None)
    parser.add_argument('--method', help='Majority call or random sampling?',
                        choices=['majority', 'random'], required=True)
    parser.add_argument('--strand-check', help='How to check for damage?',
                        choices=['USER', 'non-USER_term3', 'non-USER_all',
                        'none'], default='none')

    # if there were no arguments supplied to the function, use sys.argv
    # (skipping the first element, i.e. the name of this script)
    args = parser.parse_args(sys.argv[1:] if not argv else argv)

    result = []
    with pysam.AlignmentFile(args.bam) as bam:
        ref = pysam.FastaFile(args.ref)
        bed = BedTool(args.bed)

        # iterate through BED records and perform a pileup for each of them
        for site in bed:
            ref_base = ref.fetch(site.chrom, site.start, site.end)

            pileup_at_site = get_pileup_info(site.chrom, site.start, site.end, bam, ref_base)

            if args.strand_check:
                pileup_at_site = filter_out_damage(pileup_at_site, args.strand_check)

            if len(pileup_at_site) > 0:
                called_base = call_base(pileup_at_site, args.method)
                result.append((site.chrom, site.start, site.end, ref_base, called_base))

    for chrom, start, end, ref_base, called_base in result:
        print(chrom, start, end, ref_base, called_base, sep="\t")

    return result

if __name__ == "__main__":
    main()
