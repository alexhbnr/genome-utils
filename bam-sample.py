import random
import argparse
import sys
from collections import Counter
from itertools import chain
import pysam
from pybedtools import BedTool


def check_position(pos, read_len, is_reverse, strand_check):
    '''Test if a position in a read is potentially informative of ancient
    DNA damage (C->T or G->A substitutions) given a library preparation
    method.
    '''
    if strand_check == 'USER':
        if not is_reverse and (pos == 0 or read_len - pos <= 2): return True
        if     is_reverse and (pos  < 2 or read_len - pos == 1): return True

    elif strand_check == 'non-USER_term3' and (pos < 3 or (read_len - pos <= 3)):
        return True

    elif strand_check == 'non-USER_all':
        return True

    else:
        return False


def damage_at_site(pileup_info, strand_check):
    '''Ignore the base in the pileup of reads in case that:
       A) reference is C
          and:
           - USER treated: read has T on the first position or on the last two
           - non-USER treated:
               a) read has T on first three or last three positions
               b) read has T anywhere on the forward read
       B) reference is G
          and:
           - USER treated: read has A on the first two positions or on the last
           - non-USER treated:
               a) read has A at first three or last three positions
               b) read has A anywhere on the reverse read
    '''
    ref_base, read_base, pos_in_read, read_len, reverse_strand = pileup_info

    # if there is a C->T (on forward strand) or G->A (on reverse strand)
    # substitution at this site...
    if ((ref_base == 'C' and read_base == 'T' and not reverse_strand) or \
        (ref_base == 'G' and read_base == 'A' and     reverse_strand)):
        # ... check if it occured on a position likely to carry aDNA damage
        return check_position(pos_in_read, read_len, reverse_strand,
                              strand_check)
    return False

 
def filter_out_damage(pileup_column, strand_check):
    '''Filter out bases in a given pileup column that are likely result
    of DNA damage (C->T on forward strand, G->A on reverse strand).
    '''
    return [pileup_info for pileup_info in pileup_column
                        if not damage_at_site(pileup_info, strand_check)]


def call_base(pileup_info, sampling_method):
    """Return the most frequently occuring element of a list."""
    bases = [base for _, base, _, _, _ in pileup_info]

    if sampling_method == 'majority':
        counts = Counter(bases).most_common()
        # take all bases with the highest count
        max_freq = max(c[1] for c in counts)
        bases = [c[0] for c in counts if c[1] == max_freq]

    return random.choice(bases)


def bases_in_column(column, ref_base):
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
            pileup.append((ref_base,
                           read_base,
                           pos_in_read,
                           read_len,
                           is_reverse))

    return pileup


def sample_bases(bam, ref, sampling_method, strand_check=None, chrom=None,
                 start=None, end=None):
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

            pileup_bases = bases_in_column(col, ref_base)

            if strand_check:
                pileup_bases = filter_out_damage(pileup_bases, strand_check)

            # if there is any base in the pileup left, call one allele
            if len(pileup_bases) > 0:
                called_base = call_base(pileup_bases, sampling_method)
                sampled_bases.append((col.reference_name,
                                      col.pos + 1,
                                      ref_base, called_base))

    return sampled_bases


def sample_in_regions(bam, bed, ref, sampling_method, strand_check=None):
    '''Sample alleles from the BAM file at each position specified in a BED
    file. Return the result as a list of tuples in the form of
    (chromosome, position, ref_base, called_base).
    '''
    sampled_bases = []

    for region in bed:
        called_bases = sample_bases(bam, ref, sampling_method, strand_check,
                                    region.chrom, region.start, region.end)
        sampled_bases.append(called_bases)

    return chain.from_iterable(sampled_bases)


def print_vcf(results, sample_name, handle):
    '''Print the results in a VCF format.'''
    print('##fileformat=VCFv4.1\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}'.
          format(sample=sample_name), file=handle)

    row = '{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}'
    for chrom, pos, ref, called in results:
        alt, gt = ('.', 0) if ref == called else (called, 1)
        print(row.format(chrom=chrom, pos=pos, ref=ref, alt=alt, gt=gt),
              file=handle)


def print_bed(results, handle):
    '''Print the results in a BED format.'''
    row = '{chrom}\t{start}\t{end}\t{allele}'
    for chrom, pos, ref, called in results:
        print(row.format(chrom=chrom,
                         start=pos -1,
                         end=pos,
                         allele=called), file=handle)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Sample alleles from a given'
        ' BAM file based on pileup of reads, either by drawing random bases or'
        ' by performing a majority call at each position (from the whole BAM'
        ' or limited to regions specified by a BED file).')
    parser.add_argument('--bam', help='BAM file to sample from', required=True)
    parser.add_argument('--chrom', help='Chromosome to sample')
    parser.add_argument('--bed', help='BED file with coordinates of regions'
                        '/sites to sample')
    parser.add_argument('--ref', help='FASTA reference', required=True)
    parser.add_argument('--output', help='Name of the output file '
                        '(direct output to stdout if missing)', default=None)
    parser.add_argument('--format', help='Output as VCF or BED?',
                        choices=['VCF', 'BED'], required=True)
    parser.add_argument('--sample-name', help='Sample name to put in VCF')
    parser.add_argument('--method', help='How to sample alleles?',
                        choices=['majority', 'random'], required=True)
    parser.add_argument('--strand-check', help='How and where to check for '
                        'damage? If not specified, no checks are performed.',
                        choices=['USER', 'non-USER_term3', 'non-USER_all'],
                        default=None)

    # if there were no arguments supplied to the main function, use sys.argv
    # (skipping the first element, i.e. the name of this script)
    args = parser.parse_args(argv if argv else sys.argv[1:])

    if args.format == 'VCF' and not args.sample_name:
        parser.error('Sample has to be specified when outputting to VCF')

    bam = pysam.AlignmentFile(args.bam)
    ref = pysam.FastaFile(args.ref)

    # if user specified a BED file, perform pileup on each region in that file
    if args.bed:
        bed = BedTool(args.bed)
        results = sample_in_regions(bam, bed, ref, args.method,
                                    args.strand_check)
    else: # otherwise scan the whole BAM file directly
        results = sample_bases(bam, ref, args.method, args.strand_check,
                               args.chrom)

    # output the results as specified by user
    handle = open(args.output, 'w') if args.output else sys.stdout

    if args.format == 'VCF':
        print_vcf(results, args.sample_name, handle)
    elif args.format == 'BED':
        print_bed(results, handle)

    if args.output:
        handle.close()


if __name__ == "__main__":
    main()
