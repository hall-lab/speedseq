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
bamcleanheader.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: remove illegal and malformed fields from a BAM file's header")
    parser.add_argument('-S', '--is_sam', required=False, action='store_true', help='input is SAM')
    # parser.add_argument('-H', required=False, action='store_true', help='output header and quit')
    parser.add_argument('input', nargs='?', type=str, default=None,
                        help='SAM/BAM file to inject header lines into. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    # send back the user input
    return args

# extract read group information from header of original bam
def get_clean_header(bam):
    clean_header_list = list()

    for line in bam.text.split('\n'):
        if len(line.rstrip()) == 0:
            continue
        v = line.rstrip().split('\t')

        if v[0] == "@HD":
            legal_fields = ['VN',
                            'SO']

        elif v[0] == "@SQ":
            legal_fields = ['SN',
                            'LN',
                            'AS',
                            'M5',
                            'SP',
                            'UR']

        elif v[0] == "@RG":
            legal_fields = ['ID',
                            'CN',
                            'DS',
                            'DT',
                            'FO',
                            'KS',
                            'LB',
                            'PG',
                            'PI',
                            'PL',
                            'PU',
                            'SM']

        elif v[0] == "@PG":
            legal_fields = ['ID',
                            'PN',
                            'CL',
                            'PP',
                            'DS',
                            'VN']
        elif v[0] == "@CO":
            # all fields are legal
            clean_header_list.append(line.rstrip())
            continue

        # illegal tag
        else:
            continue

        tag = dict(x.split(':',1) for x in v[1:])
        tag_clean = dict()            

        for field in tag:
            if (field in legal_fields
                and tag[field].strip() != ""):
                tag_clean[field] = tag[field]

        clean_header_list.append(v[0] + '\t'
                                 + '\t'.join(field + ":" + tag[field] for field in tag))
    return '\n'.join(clean_header_list)

# add read group info to header of new sam file
def bam_clean(bam, is_sam, header_only):
    if is_sam:
        in_bam = pysam.Samfile(bam, 'r')
    else:
        in_bam = pysam.Samfile(bam, 'rb')

    # out_bam = pysam.Samfile('-', 'w', template=in_bam)

    print get_clean_header(in_bam)

    if not header_only:
        for al in in_bam:
            print al

    # # this code leads to pipeing errors
    # if not header_only:
    #     for al in in_bam:
    #         out_bam.write(al)

    # out_bam.close()
    # close the input file

    in_bam.close()

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # clean the header
    bam_clean(args.input, args.is_sam, True)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
