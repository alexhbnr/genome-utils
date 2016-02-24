##  bam-sample.py

Sample alleles from a given BAM file based on pileup of reads, either by drawing random bases or by performing a majority call at each position (from the whole BAM or limited to regions specified by a BED file).

If necessary, any potential ancient DNA damage can be taken care of accordingly based on characteristics of a library preparation method:

* For USER treated libraries, sites carrying C → T substitutions (on forward reads) or G → A substitutions (on reverse reads) on the first position and on the last two positions of reads will be ignored during sampling.

* For non-USER treated libraries:

  - Option 1. Ignore C → T substitutions on first three and last three positions of forward reads. Similarly, ignore G → A substitutions on terminal three positions of reverse reads.
  
  - Option 2. Ignore C → T or G → A substitutions throughout the whole length of forward or reverse reads respectively.
