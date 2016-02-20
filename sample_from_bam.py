import random
import argparse
import sys
from collections import Counter
from itertools import chain
import pysam
from pybedtools import BedTool


def check_position(pos, read_len, reverse_strand, library_prep):
    '''Test if a position in a read is potentially informative of ancient
    DNA damage (C->T or G->A substitutions) given a library preparation
    method.
    '''
    if library_prep == 'USER':
        if not reverse_strand and ((pos == 0) or (read_len - pos <= 2)): return True
        if     reverse_strand and ((pos  < 2) or (read_len - pos == 1)): return True

    elif library_prep == 'non-USER_term3' and ((pos < 3) or (read_len - pos <= 3)):
        return True

    elif library_prep == 'non-USER_all':
        return True

    else:
        return False


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

    # if there is a potential C->T or G->A substitution, check if
    # it occured on a position that is likely to carry ancient DNA damage
    if ((ref_base == 'C' and read_base == 'T' and not reverse_strand) or \
        (ref_base == 'G' and read_base == 'A' and     reverse_strand)):
        return check_position(pos_in_read, read_len, reverse_strand, library_prep)

    return False

 
def filter_out_damage(pileup_column, library_prep):
    '''Filter out bases in a given pileup column that are likely result
    of DNA damage (C->T on forward strand, G->A on reverse strand).
    '''
    return [site for site in pileup_column if not damage_at_site(site, library_prep)]


def call_base(pileup_info, sampling_method):
    """Return the most frequently occuring element of a list."""
    bases = [base for base, _, _, _ in pileup_info]

    if sampling_method == 'majority':
        counts = Counter(bases).most_common()
        # take all bases with the highest count
        max_freq = max(c[1] for c in counts)
        bases = [c[0] for c in counts if c[1] == max_freq]

    return random.choice(bases)


def bases_in_column(column):
    '''Return a list of bases in a given pileup column.
    '''
    pileup = []

    # walk through all reads overlapping the current column
    # and accumulate bases at that position
    for pileup_read in column.pileups:
        # skip deletions
        if pileup_read.is_del: continue

        pos_in_read = pileup_read.query_position
        read_len = pileup_read.alignment.query_length
        read_base = pileup_read.alignment.query_sequence[pos_in_read]
        is_reverse = pileup_read.alignment.is_reverse

        if read_base in "ACGT": 
            pileup.append((read_base,
                           pos_in_read,
                           read_len,
                           is_reverse))

    return pileup


def sample_bases(bam, ref, sampling_method, strand_check=None, chrom=None, start=None, end=None):
    '''Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    '''
    sampled_bases = []

    for col in bam.pileup(chrom, start, end):
        # if coordinates were specified, check first if a current column
        # lies within this region (pysam pileui return whole reads overlapping
        # a requested region, not just bases in this region)
        if chrom and start and end and not (start <= col.pos < end):
            continue
        else:
            ref_base = ref.fetch(col.reference_name, col.pos, col.pos + 1)

            pileup_bases = bases_in_column(col)

            if strand_check:
                pileup_bases = filter_out_damage(pileup, strand_check)

            # if there is any base in the pileup left, call one allele
            if len(pileup_bases) > 0:
                called_base = call_base(pileup_bases, sampling_method)
                sampled_bases.append((col.reference_name,
                                      col.pos + 1,
                                      ref_base, called_base))

    return sampled_bases


def sample_in_regions(bam, bed, ref, sampling_method, library_prep=None):
    '''Sample alleles from the BAM file at each position specified in a BED
    file. Return the result as a list of tuples in the form of
    (chromosome, position, ref_base, called_base).
    '''
    sampled_bases = []

    for region in bed:
        called_bases = sample_bases(bam, ref, sampling_method, library_prep,
                                   region.chrom, region.start, region.end)
        sampled_bases.append(called_bases)

    return chain.from_iterable(sampled_bases)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Sample bases from BAM file')
    parser.add_argument('--bam', help='BAM file to sample from', required=True)
    parser.add_argument('--bed', help='BED file with coordinates of sites'
                        'to sample at')
    parser.add_argument('--ref', help='FASTA reference', required=True)
    parser.add_argument('--output', help='Name of the output file '
                        '(direct output to stdout if missing)', default=None)
    parser.add_argument('--sampling-method', help='Majority or random call?',
                        choices=['majority', 'random'], required=True)
    parser.add_argument('--strand-check', help='How to check for damage '
                        '(this is determined by library preparation method)',
                        choices=['USER', 'non-USER_term3', 'non-USER_all'],
                        default=None)

    # if there were no arguments supplied to the main function, use sys.argv
    # (skipping the first element, i.e. the name of this script)
    args = parser.parse_args(argv if argv else sys.argv[1:])

    bam = pysam.AlignmentFile(args.bam)
    ref = pysam.FastaFile(args.ref)

    # if user specified a BED file, perform pileup on each region in that file
    if args.bed:
        bed = BedTool(args.bed)
        results = sample_in_regions(bam, bed, ref, args.sampling_method, args.strand_check)
    else: # otherwise scan the whole BAM file directly
        results = sample_bases(bam, ref)

    for chrom, pos, ref_base, called_base in results:
        print(chrom, pos, ref_base, called_base, sep='\t')

    return results

if __name__ == "__main__":
    main()
