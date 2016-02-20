import random
import argparse
import sys
from collections import Counter
import pysam
from pybedtools import BedTool


def terminal_base(pos, read_len, is_forward, library_prep):
    '''Check if a given position is on the beginning or the end a read.
    What exactly is a 'beginning' and 'end' of a read is determined by
    strand orientation and the library preparation method which influences
    the patterns of substitutions (C->T, G->A) in damaged ancient DNA.
    '''
    is_terminal = False

    if library_prep == 'USER':
        if     is_forward and ((pos == 0) or (read_len - pos <= 2)): is_terminal = True
        if not is_forward and ((pos  < 2) or (read_len - pos == 1)): is_terminal = True
    if library_prep == 'non-USER_term3':
        if (pos < 3) or (read_len - pos <= 3): is_terminal = True

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
    of DNA damage (C->T on forward strand, G->A on reverse strand).
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


def call_base(pileup_info, sampling_method):
    """Return the most frequently occuring element of a list."""
    bases = [base for _, base, _, _, _ in pileup_info]

    if sampling_method == 'majority':
        counts = Counter(bases).most_common()
        # take all bases with the highest count
        max_freq = max(c[1] for c in counts)
        bases = [c[0] for c in counts if c[1] == max_freq]

    return random.choice(bases)


def scan_bam(bam_file, bed_file, ref_file, library_prep, sampling_method):
    '''Sample alleles from the BAM file at each position specified in a BED
    file. Return the result as a list of tuples in the form of
    (chromosome, position, ref_base, called_base).
    '''
    result = []

    with pysam.AlignmentFile(bam_file) as bam:
        ref = pysam.FastaFile(ref_file)
        bed = BedTool(bed_file)

        # iterate through BED records and perform a pileup for each of them
        for site in bed:
            ref_base = ref.fetch(site.chrom, site.start, site.end)

            pileup_at_site = get_pileup_info(site.chrom, site.start, site.end, bam, ref_base)

            if library_prep:
                pileup_at_site = filter_out_damage(pileup_at_site, library_prep)

            if len(pileup_at_site) > 0:
                called_base = call_base(pileup_at_site, sampling_method)
                result.append((site.chrom, site.end, ref_base, called_base))

    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description='Sample bases from BAM file')
    parser.add_argument('--bam', help='BAM file to sample from', required=True)
    parser.add_argument('--bed', help='BED file with coordinates of sites'
                        'to sample at', required=True)
    parser.add_argument('--ref', help='FASTA reference', required=True)
    parser.add_argument('--output', help='Name of the output file '
                        '(direct output to stdout if missing)', default=None)
    parser.add_argument('--sampling-method', help='Majority or random call?',
                        choices=['majority', 'random'], required=True)
    parser.add_argument('--strand-check', help='How to check for damage?',
                        choices=['USER', 'non-USER_term3', 'non-USER_all',
                        'none'], default='none')

    # if there were no arguments supplied to the main function, use sys.argv
    # (skipping the first element, i.e. the name of this script)
    args = parser.parse_args(argv if argv else sys.argv[1:])

    # list to accumulate tuples of (chrom, pos, ref_base, called_base)
    result = scan_bam(args._bam, args._bed, args._ref, args._library_prep,
                      args._sampling_method)

    for chrom, pos, ref_base, called_base in result:
        print(chrom, pos, ref_base, called_base, sep="\t")

    return result

if __name__ == "__main__":
    main()
