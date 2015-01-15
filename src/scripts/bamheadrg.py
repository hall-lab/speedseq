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
bamheadrg.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Inject readgroup info")
    parser.add_argument('-r', '--readgroup', default=None, required=False, help='Read group(s) to extract (comma separated)')
    parser.add_argument('-d', '--donor', type=str, required=True, help='Donor BAM/SAM file to extract read group info')
    parser.add_argument('-S', '--donor_is_sam', required=False, action='store_true', help='Donor file is SAM')
    parser.add_argument('recipient', nargs='?', type=argparse.FileType('r'), default=None,
                        help='SAM file to inject header lines into. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.recipient == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.recipient = sys.stdin

    # send back the user input
    return args

# extract read group information from header of original bam
def extract_rg_info(donor, donor_is_sam, rgs_to_extract):
    if donor_is_sam:
        bam = pysam.Samfile(donor, 'r')
    else:
        bam = pysam.Samfile(donor, 'rb')
    rg_out = list()
    for readgroup in bam.header['RG']:
        if not rgs_to_extract or readgroup['ID'] in rgs_to_extract:
            rg_out.append(readgroup)
    bam.close()
    return rg_out

# add read group info to header of new sam file
def bamheadrg(recipient, rg_out):
    in_header = True
    for line in recipient:
        if in_header:
            if line[0] != '@':
                for readgroup in rg_out:
                    print '@RG\t' + '\t'.join([':'.join((t,readgroup[t])) for t in readgroup])
                in_header = False
        print line.rstrip()
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    if args.readgroup:
        rgs_to_extract = args.readgroup.split(',')
    else:
        rgs_to_extract = None

    # extract specified readgroups from original bam file
    rg_out = extract_rg_info(args.donor, args.donor_is_sam, rgs_to_extract)

    # add extracted readgroups to new sam file
    bamheadrg(args.recipient, rg_out)

    # close the input file
    args.recipient.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
