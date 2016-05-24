#!/usr/bin/env python
import sys, os, os.path, subprocess, time, re
import argparse
from argparse import RawTextHelpFormatter

__author__ = "Allison Regier (aregier@wustl.edu)"
__version__ = "$Revision: 0.0.2 $"
__date__ = "$Date: 2016-05-09 14:53 $"

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
cnvnator\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: SpeedSeq wrapper for CNVnator v0.3.2")
    parser.add_argument('-w', '--window', type=str, required=False, default='100', help='window size in base pairs [100]')
    parser.add_argument('-b', '--bam', required=True, help='input bam file name')
    parser.add_argument('-o', '--output', required=True, help='output variant file name')
    parser.add_argument('-c', '--chroms', required=True, help='path to chromosome files')
    parser.add_argument('-g', '--genome', required=False, default='GRCh37', help='genome build [GRCh37]')
    parser.add_argument('--cnvnator', required=False, default='cnvnator-multi', help='path to cnvnator-multi binary')
    parser.add_argument('-T', '--tempdir', type=str, required=False, default='temp', help='temp directiory [./temp]')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args

# get root file name
def get_root_fn(bam_fn):
	return TEMPDIR + '/' + os.path.split(bam_fn)[1] + ROOT_EXT
# end of get root file name


# get hist file name
def get_hist_fn(bam_fn):
	return TEMPDIR + '/' + os.path.split(bam_fn)[1] + HIST_EXT + ROOT_EXT
# end of get hist file name


# get list of chromosomes from the BAM file header
def get_chroms_list(bam_fn):
	proc = subprocess.Popen(['samtools', 'view', '-H', bam_fn], stdout = subprocess.PIPE)
	(dout, derr) = proc.communicate()	
	chroms_list = []
	lines = dout.split('\n')
	for line in lines:
		pieces = line.split()
		if len(pieces) != 3: continue
		if pieces[0] == "@SQ":
			chrm = pieces[1][pieces[1].find(":")+1:]
			chroms_list.append(chrm)	
	return chroms_list
# end of chromosomes list


# run_cnvnator
def run_partition(bin_size, root_fn, chroms):
	devnull = open(os.devnull, 'w')
	ret = subprocess.call([CNVNATOR, '-root', root_fn, '-partition', bin_size], stdout = devnull)
	devnull.close()
	return ret
# end of run_cnvnator


# run tree, hist, and stats
def run_hist_stats(bin_size, bam_fn, chroms_dir):
	root_fn = get_root_fn(bam_fn)
	hist_fn = get_hist_fn(bam_fn)

	print "===== Running histograms on input data for input bin size"
	ret = subprocess.call([CNVNATOR, '-his', bin_size, '-d', chroms_dir, '-root', root_fn, '-outroot', hist_fn]) 
	if ret != 0:
		print "Error computing histograms (input bin size)."
		return ret
	print "===== Running stats on input data for input bin size"
	ret = subprocess.call([CNVNATOR, '-stat', bin_size, '-root', hist_fn]) 
	if ret != 0:
		print "Error computing histograms (input bin size)."
		return ret
        # no need to duplicate hist and stat if the input bin size was 1000
        if bin_size == "1000":
            return 0
        
        #print "===== Running histograms on input data for bin size 1000"
        ret = subprocess.call([CNVNATOR, '-his', '1000', '-d', chroms_dir, '-root', root_fn, '-outroot', hist_fn]) 
        if ret != 0:
            print "Error computing histograms (bin size 1000)."
            return ret
        #print "===== Running stats on input data for bin size 1000"
        ret = subprocess.call([CNVNATOR, '-stat', '1000', '-root', hist_fn]) 
        if ret != 0:
            print "Error computing stats (bin size 1000)."
            return ret
        return 0
# end of run tree, hist, stats


# run calls
def run_calls(bin_size, hist_fn, out_fn):
	print "===== Running calls on input data"
	f = open(out_fn + '.txt', 'w')
	ret = subprocess.call([CNVNATOR, '-call', bin_size, '-root', hist_fn], stdout = f) 
	f.close()
	if ret != 0:
		print "Error computing calls."
	return ret
# end of run calls


# run tree
def run_tree(bam_fn, genome, chroms):
	print "===== Running tree on input data"
	#for x in chroms:
	#	ret = subprocess.call([CNVNATOR, '-root', get_root_fn(bam_fn), '-genome', genome, '-tree', bam_fn, '-unique', '-chrom', x])
	#	if ret != 0:
	#		print "Error in tree creation chromosome x"
	#		return ret
    #return 0
	print "Running on bam %s" % bam_fn
	sep=" "
	ret = subprocess.call([CNVNATOR, '-root', get_root_fn(bam_fn), '-genome', genome, '-tree', bam_fn, '-unique', '-chrom', sep.join(chroms)])
	return ret
# end of run tree


# make a bedgraph file
def mk_graph_file(out_fn):
	bedgraph_fn = out_fn + ".bed"
	f = open(out_fn + '.txt', 'r')
	fdata = f.readlines()
	f.close()
	
	nf = []
	prev_chr = ""
	prev_end = 0
	for x in fdata:
  		pieces = x.split()
  		idx = pieces[1].find(':')
  		idx2 = pieces[1].find('-')
  		chr = pieces[1][0:idx]
#		# remove any 'chr' prefixes from the chrom name
#		if chr.startswith('chr'): chr = chr[3:]
  		start = int(pieces[1][idx+1:idx2])
  		end = int(pieces[1][idx2+1:])
  		prev_chr, prev_end = chr, end
  		line = "%s\t%d\t%d\t%s" % (chr, start, end, pieces[3])
  		nf.append(line)
	nfstr = '\n'.join(nf)
	
	f = open(bedgraph_fn, 'w')
	f.write(nfstr)
	f.close()
# end of make bedgraph file


# main
if __name__ == "__main__":
	args = get_args()

        CNVNATOR = args.cnvnator
        TEMPDIR = args.tempdir
        ROOT_EXT = ".root"
        HIST_EXT = ".hist"
	
        # create temp directory if it doesn't exist
        try:
            os.stat(TEMPDIR)
        except:
            os.mkdir(TEMPDIR)

        # build chroms_list
	chroms_list = get_chroms_list(args.bam)
	if len(chroms_list) == 0:
		print "No chromosomes found in BAM file."
		sys.exit(1)
	print "Processing data from the following chromosomes: %s" % str(chroms_list)
	
	# run tree
	if run_tree(args.bam, args.genome, chroms_list) != 0:
		sys.exit(1)
        
	# run hist and stats
	if run_hist_stats(args.window, args.bam, args.chroms) != 0:
		sys.exit(1)

	# # run partition
	hist_fn = get_hist_fn(args.bam)	
	if run_partition(args.window, hist_fn, chroms_list) != 0:
		sys.exit(1)

	# run calls
        print args.output
	if run_calls(args.window, hist_fn, args.output) != 0:
		sys.exit(1)
	mk_graph_file(args.output)
	sys.exit(0)
# end of main

