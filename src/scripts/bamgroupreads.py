#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
import string
from string import *

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-12-15 11:43 $"

def bamgroupreads(bamfile, readgroup, reset_dups, fix_flags, is_sam, bam_out, uncompressed_out):
    # set input file
    if bamfile == None: 
        if is_sam:
            in_bam = pysam.Samfile("-", "r")
        else:
            in_bam = pysam.Samfile('-', 'rb')
    else:
        if is_sam:
            in_bam = pysam.Samfile(bamfile, 'r')
        else:
            in_bam = pysam.Samfile(bamfile, "rb")

    # set output file
    if uncompressed_out:
        out_bam = pysam.Samfile('-', 'wbu', template=in_bam)
    elif bam_out:
        out_bam = pysam.Samfile('-', 'wb', template=in_bam)
    else:
        out_bam = pysam.Samfile('-', 'wh', template=in_bam)
        

    # parse readgroup string
    try:
        rg_list = readgroup.split(',')
    except AttributeError:
        rg_list = None

    d = {}
    for al in in_bam.fetch():
        # must be in a user specified readgroup
        if rg_list and al.opt('RG') not in rg_list:
            continue

        # add read name to dictionary if not already there
        key = al.qname
        if key not in d:
            d.setdefault(key,Namegroup(al))
        # print matched read pairs
        else:
            d[key].add_alignment(al)
            if d[key].is_complete():
                for al in d[key].alignments:
                    if reset_dups:
                        # unset the duplicate flag
                        al.is_duplicate = 0
                    if fix_flags:
                        # fix the secondary mate flag
                        proper_pair = False
                        duplicate = False
                        read1_unmapped = False
                        read2_unmapped = False

                        # gather info on the read cluster flags
                        for flagcheck in d[key].alignments:
                            if flagcheck.is_proper_pair:
                                proper_pair = True
                            if flagcheck.is_duplicate:
                                duplicate = True
                            if flagcheck.is_secondary:
                                continue
                            if flagcheck.is_read1:
                                read1_unmapped = flagcheck.is_unmapped
                            elif flagcheck.is_read2:
                                read2_unmapped = flagcheck.is_unmapped

                        # set new info on the read cluster
                        if al.is_read1:
                            al.mate_is_unmapped = read2_unmapped
                        elif al.is_read2:
                            al.mate_is_unmapped = read1_unmapped
                        al.is_proper_pair = proper_pair
                        al.is_duplicate = duplicate
                    out_bam.write(al)
                del d[key]
    if len(d) != 0:
        sys.stderr.write('Error: %s unmatched name groups\n' % len(d))
        exit(1)

# ============================================
# functions
# ============================================

# class that holds reads from a sequence fragment
class Namegroup():
    def __init__(self, al):
        self.alignments = list()
        self.name = al.qname
        self.sa = 0
        self.num_prim = 0
        self.add_alignment(al)

    def add_alignment(self, al):
        self.alignments.append(al)
        if not al.is_secondary:
            self.num_prim += 1
            try:
                self.sa += len(al.opt('SA').rstrip(';').split(';'))
                # print self.sa
            except KeyError:
                pass

    def is_complete(self):
        return self.num_prim == 2 and len(self.alignments) == self.sa + 2

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Group BAM file by read IDs without sorting")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-r', '--readgroup', metavar='STR', default=None, required=False, help='Read group(s) to extract (comma separated)')
    parser.add_argument('-d', '--reset_dups', required=False, action='store_true', help='Reset duplicate flags')
    parser.add_argument('-f', '--fix_flags', required=False, action='store_true', help='Fix mate flags for secondary reads')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-b', required=False, action='store_true', help='Output BAM format')
    parser.add_argument('-u', required=False, action='store_true', help='Output uncompressed BAM format (implies -b)')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args

# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    args = get_args()
    bamgroupreads(args.input, args.readgroup, args.reset_dups, args.fix_flags, args.S, args.b, args.u)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
    
