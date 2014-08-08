#!/usr/bin/env python

import pysam
import argparse, sys
import math, time, re
import multiprocessing
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-04-28 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
svgt\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Compute genotype of structural variants based on breakpoint depth")
    parser.add_argument('-B', '--bam', type=str, required=True, help='BAM file(s), comma-separated if genotyping multiple BAMs')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args


def get_mate(bam, read):
    if read.is_proper_pair and read.tid == read.rnext and not read.mate_is.unmapped:
        # chrom = bam.getrname(read.tid)
        # for mate in bam.fetch(chrom, read.pnext, read.pnext + 1):
        #     if mate.qname == read.qname:
        #         break
        # return mate
        pointer = bam.tell()

        try:
            mate = bam.mate(read)
        finally:
            bam.seek(pointer)

        # print read.pnext, read.
        return mate

def default_get_mate(bam, read):
    return bam.mate(read)

# primary function
def sv_genotype(bamstr):
    bam = pysam.Samfile(bamstr, 'rb')

    lim = 5
    counter = 0
    # for read in bam.fetch('8', 43099856, 43100401):
    # for read in bam.fetch('8', 43099856, 43099857):
    for read in bam.fetch('8', 43454401, 43454502):
        if counter <= lim:
            print read
            # print default_get_mate(bam, read)
            print get_mate(bam, read)
            counter += 1
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    sv_genotype(args.bam)


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
