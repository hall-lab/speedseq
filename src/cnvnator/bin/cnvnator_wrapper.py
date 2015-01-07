#!/usr/bin/env python
import sys, os, os.path, subprocess, time 
import argparse
from argparse import RawTextHelpFormatter

__author__ = "David Rose (dbr3d@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-08-14 14:53 $"

before_part_cp = """
int before_partition_copy(TString fname_in  = "his.HG00702.root",
			  TString fname_out = "22.root",
			  TString bin       = "1000",
			  TString chr       = "22")
{
  TString dirname = "bin_"; dirname += bin;

  TFile      *file_in  = new TFile(fname_in, "READ");
  if (!file_in) {
    cerr<<"No file '"<<file_in<<"' found."<<endl;
    return 1;
  }
  TDirectory *dir_in   = file_in->GetDirectory(dirname);
  if (!dir_in) {
    cerr<<"No dir '"<<dirname<<"' in the input file."<<endl;
    return 1;
  }

  TFile      *file_out = new TFile(fname_out,"UPDATE");
  TDirectory *dir_out  = file_out->GetDirectory(dirname);
  if (!dir_out) {
    cerr<<"No dir '"<<dirname<<"' in the output file."<<endl;
    cerr<<"Creating ..."<<endl;
    dir_out = file_out->mkdir(dirname);
  }
  dir_out->cd();

  TString prefix = "his_rd_p_"; prefix += chr; prefix += "_"; prefix += bin;
  TH1 *his = NULL;

  TString hisname = prefix + "_GC";
  dir_in->GetObject(hisname,his);
  his->Write(his->GetName(),TObject::kOverwrite);

  hisname = "rd_p_GC_" + bin;
  dir_in->GetObject(hisname,his);
  his->Write(his->GetName(),TObject::kOverwrite);

  hisname = "rd_p_xy_GC_" + bin;
  dir_in->GetObject(hisname,his);
  his->Write(his->GetName(),TObject::kOverwrite);

  file_out->Close();
  file_in->Close();
  return 0;
}
"""

after_part_cp = """
int after_partition_copy(TString fname_in  = "22.root",
			 TString fname_out = "his.22.root", // here put name of master file
			 TString bin       = "1000",
			 TString chr       = "22")
{
  TString dirname = "bin_"; dirname += bin;

  TFile      *file_in  = new TFile(fname_in, "READ");
  if (!file_in) {
    cerr<<"No file '"<<file_in<<"' found."<<endl;
    return 1;
  }
  TDirectory *dir_in   = file_in->GetDirectory(dirname);
  if (!dir_in) {
    cerr<<"No dir '"<<dirname<<"' in the input file."<<endl;
    return 1;
  }

  TFile      *file_out = new TFile(fname_out,"UPDATE");
  TDirectory *dir_out  = file_out->GetDirectory(dirname);
  if (!dir_out) {
    cerr<<"No dir '"<<dirname<<"' in the output file."<<endl;
    cerr<<"Creating ..."<<endl;
    dir_out = file_out->mkdir(dirname);
  }
  dir_out->cd();

  TString prefix = "his_rd_p_"; prefix += chr; prefix += "_"; prefix += bin;
  TH1 *his = NULL;

  TString hisname = prefix + "_GC_l1";
  dir_in->GetObject(hisname,his);
  if (his) his->Write(his->GetName(),TObject::kOverwrite);
  his = NULL;

  hisname = prefix + "_GC_l2";
  dir_in->GetObject(hisname,his);
  if (his) his->Write(his->GetName(),TObject::kOverwrite);
  his = NULL;

  hisname = prefix + "_GC_l3";
  dir_in->GetObject(hisname,his);
  if (his) his->Write(his->GetName(),TObject::kOverwrite);
  his = NULL;

  hisname = prefix + "_partition_GC";
  dir_in->GetObject(hisname,his);
  if (his) his->Write(his->GetName(),TObject::kOverwrite);
  his = NULL;
  
  file_out->Close();
  file_in->Close();

  return 0;
}
"""


def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
cnvnator\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: SpeedSeq parallelized implementation of CNVnator v0.3 (Gerstein lab)")
    parser.add_argument('-t', '--threads', type=int, required=False, default=1, help='number of threads')
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
	ret = proc.poll()
	while ret is None:
		time.sleep(2)
		ret = proc.poll()
	(dout, derr) = proc.communicate()	
	if ret != 0:
		print "Error viewing / processing header from specified BAM input file: %s" % derr
		return []
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
	# write partition copy scripts if they don't already exist
	if not os.path.exists(BEFORE_PART_FN):	
		with open(BEFORE_PART_FN, 'w') as f: f.write(before_part_cp)
	if not os.path.exists(AFTER_PART_FN):
		with open(AFTER_PART_FN, 'w') as f: f.write(after_part_cp)
	
	# before_parition_copy and after_paritition_copy are called in two ways, depending on chromosome
	def get_part_cp_chrom(chrom):
                return chrom
	# end

	# copy data from the ROOT file to a standalone data file for each chromosome
	devnull = open(os.devnull, 'w')
	for x in chroms:
		print "Extracting input data for chrom %s..." % x
		proc = subprocess.call(['root', '-b', '-q', '%s(\"%s\",\"%s/%s.root\", \"%s\", \"%s\")' % (BEFORE_PART_FN, root_fn, TEMPDIR, x, bin_size, get_part_cp_chrom(x))], stdout = devnull, stderr = devnull)
		if proc != 0:
			print "Error: Data extraction (before_partition_copy) failed for chrom %s" % x
                        exit(1)
	
	# spawn off up to max_procs copies of CNVNATOR to partition each chromosome
	part_dict = {}
	proc_idx = 0
	while proc_idx < len(chroms): 
		x = chroms[proc_idx]
		part_dict[x] = subprocess.Popen([CNVNATOR, '-root', '%s/%s.root' % (TEMPDIR, x), '-partition', bin_size, '-chrom', x], stdout = devnull)
		print ("Now partitioning chrom %s" % x)
		proc_idx += 1
		if proc_idx >= MAX_PROCS: break

	# check on the CNVNATOR procs
	# this could take a while
	ret = 0
	while True:
		if len(part_dict) == 0: break
		time.sleep(10)
		for x in part_dict.keys():
			val = part_dict[x].poll()
			if val == None:
				continue
			elif val == 0:
				# if we get a good return, merge the chrom data back into the root file and clean up the standalone data
				print "Partitioning succeeded for chrom %s" % x
				proc = subprocess.call(['root', '-b', '-q', '%s(\"%s/%s.root\", \"%s\", \"%s\", \"%s\")' % (AFTER_PART_FN, TEMPDIR, x, root_fn, bin_size, get_part_cp_chrom(x))], stdout = devnull)			
				os.unlink('%s/%s.root' % (TEMPDIR, x))
			else:
				print "Partitioning failed for chrom %s" % x
				ret = 1
			part_dict.pop(x)
			# since we spawn up to max_procs number of procs at once, we might have more things to run...
			if proc_idx < len(chroms):
				x = chroms[proc_idx]
				part_dict[x] = subprocess.Popen([CNVNATOR, '-root', '%s/%s.root' % (TEMPDIR, x), '-partition', bin_size, '-chrom', x], stdout = devnull)
				print ("Now partitioning chrom %s" % x)
				proc_idx += 1
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
def run_tree(bam_fn, genome):
	print "===== Running tree on input data"
	ret = subprocess.call([CNVNATOR, '-root', get_root_fn(bam_fn), '-genome', genome, '-ptree', bam_fn, '-unique', '-nthreads', str(MAX_PROCS)])
	if ret != 0:
		print "Error in tree creation"
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
		# remove any 'chr' prefixes from the chrom name
		if chr.startswith('chr'): chr = chr[3:]
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
        BEFORE_PART_FN = TEMPDIR + "/" + "before_partition_copy.cpp"
        AFTER_PART_FN = TEMPDIR + "/" + "after_partition_copy.cpp"
        MAX_MAX_PROCS = 48
        MAX_PROCS = 32
        ROOT_EXT = ".root"
        HIST_EXT = ".hist"
	
        # create temp directory if it doesn't exist
        try:
            os.stat(TEMPDIR)
        except:
            os.mkdir(TEMPDIR)

        # ensure threads doesn't exceed MAX_MAX_PROCS
	try:
		MAX_PROCS = args.threads
		if MAX_PROCS <= 0 or MAX_PROCS > MAX_MAX_PROCS:
			raise Exception()	
	except:
		print "Number of threads must be a number between 1 and %d." % MAX_MAX_PROCS
		sys.exit(1)
	
        # build chroms_list
	chroms_list = get_chroms_list(args.bam)
	if len(chroms_list) == 0:
		print "No chromosomes found in BAM file."
		sys.exit(1)
	print "Processing data from the following chromosomes: %s" % str(chroms_list)
	
	# run tree
	if run_tree(args.bam, args.genome) != 0:
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

