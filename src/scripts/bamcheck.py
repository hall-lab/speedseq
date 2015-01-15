#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import pysam

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamcheck.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: check BAM for EOF")
    parser.add_argument('bam', type=str, help='BAM file to check for EOF')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # try:
    pysam.Samfile(args.bam, "rb")
    # except:
    #     [W::bam_hdr_read]


    # close the input file
    # args.bam.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
