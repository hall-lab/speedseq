#!/usr/bin/env python
import sys, os, os.path, subprocess, time 

before_part_cp = """
int before_partition_copy(TString fname_in  = "his.HG00702.root",
			  TString fname_out = "chr22.root",
			  TString bin       = "1000",
			  TString chr       = "chr22")
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
int after_partition_copy(TString fname_in  = "chr22.root",
			 TString fname_out = "his.chr22.root", // here put name of master file
			 TString bin       = "1000",
			 TString chr       = "chr22")
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

BEFORE_PART_FN = "before_partition_copy.cpp"
AFTER_PART_FN = "after_partition_copy.cpp"
MAX_PROCS = 128
ROOT_EXT = ".root"
HIST_EXT = ".hist"


# get root file name
def get_root_fn(bam_fn):
	return os.path.split(bam_fn)[1] + ROOT_EXT
# end of get root file name


# get hist file name
def get_hist_fn(bam_fn):
	return os.path.split(bam_fn)[1] + HIST_EXT + ROOT_EXT
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
		if chrom.startswith("GL"): return chrom
		return "chr%s" % chrom
	# end


	# copy data from the ROOT file to a standalone data file for each chromosome
	devnull = open(os.devnull, 'w')
	for x in chroms:
		print "Extracting input data for chrom %s..." % x
		proc = subprocess.call(['root', '-b', '-q', '%s(\"%s\",\"chr%s.root\", \"%s\", \"%s\")' % (BEFORE_PART_FN, root_fn, x, bin_size, get_part_cp_chrom(x))], stdout = devnull, stderr = devnull)
		if proc != 0:
			print "Error: Data extraction (before_partition_copy) failed for chrom %s" % x
	
	# spawn off up to max_procs copies of CNVNATOR to partition each chromosome
	part_dict = {}
	proc_idx = 0
	while proc_idx < len(chroms): 
		x = chroms[proc_idx]
		part_dict[x] = subprocess.Popen(['cnvnator', '-root', 'chr%s.root'%x, '-partition', bin_size, '-chrom', x], stdout = devnull)
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
				proc = subprocess.call(['root', '-b', '-q', '%s(\"chr%s.root\", \"%s\", \"%s\", \"%s\")' % (AFTER_PART_FN, x, root_fn, bin_size, get_part_cp_chrom(x))], stdout = devnull)			
				os.unlink('chr%s.root'%x)
			else:
				print "Partitioning failed for chrom %s" % x
				ret = 1
			part_dict.pop(x)
			# since we spawn up to max_procs number of procs at once, we might have more things to run...
			if proc_idx < len(chroms):
				x = chroms[proc_idx]
				part_dict[x] = subprocess.Popen(['cnvnator', '-root', 'chr%s.root'%x, '-partition', bin_size, '-chrom', x], stdout = devnull)
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
	ret = subprocess.call(['cnvnator', '-his', bin_size, '-d', chroms_dir, '-root', root_fn, '-outroot', hist_fn]) 
	if ret != 0:
		print "Error computing histograms (input bin size)."
		return ret
	print "===== Running stats on input data for input bin size"
	ret = subprocess.call(['cnvnator', '-stat', bin_size, '-root', hist_fn]) 
	if ret != 0:
		print "Error computing histograms (input bin size)."
		return ret
	# no need to duplicate hist and stat if the input bin size was 1000
	if bin_size == "1000":
		return 0
	
	#print "===== Running histograms on input data for bin size 1000"
	ret = subprocess.call(['cnvnator', '-his', '1000', '-d', chroms_dir, '-root', root_fn, '-outroot', hist_fn]) 
	if ret != 0:
		print "Error computing histograms (bin size 1000)."
		return ret
	#print "===== Running stats on input data for bin size 1000"
	ret = subprocess.call(['cnvnator', '-stat', '1000', '-root', hist_fn]) 
	if ret != 0:
		print "Error computing stats (bin size 1000)."
		return ret
	return 0
# end of run tree, hist, stats


# run calls
def run_calls(bin_size, hist_fn, out_fn):
	print "===== Running calls on input data"
	f = open(out_fn, 'w')
	ret = subprocess.call(['cnvnator', '-call', bin_size, '-root', hist_fn], stdout = f) 
	f.close()
	if ret != 0:
		print "Error computing calls."
	return ret
# end of run calls


# run tree
def run_tree(bam_fn, genome):
	print "===== Running tree on input data"
	ret = subprocess.call(['cnvnator', '-root', get_root_fn(bam_fn), '-genome', genome, '-tree', bam_fn, '-unique'])
	if ret != 0:
		print "Error in tree creation"
	return ret
# end of run tree


# main
if __name__ == "__main__":
	# usage: ./cnvnator.py {window size} {BAM file} {output variant file name} {path to chromosome files}
	if len(sys.argv) != 6:
		print "Usage: ./cnvator.py {window size} {BAM file name} {output variant file name} {path to chromosome files} {genome}"
		sys.exit(1)
	chroms_list = get_chroms_list(sys.argv[2])
	if len(chroms_list) == 0:
		print "No chromosomes found in BAM file."
		sys.exit(1)
	print "Processing data from the following chromosomes: %s" % str(chroms_list)

	# run tree
	if run_tree(sys.argv[2], sys.argv[5]) != 0:
		sys.exit(1)

	# run hist and stats
	if run_hist_stats(sys.argv[1], sys.argv[2], sys.argv[4]) != 0:
		sys.exit(1)

	# run partition
	hist_fn = get_hist_fn(sys.argv[2])	
	if run_partition(sys.argv[1], hist_fn, chroms_list) != 0:
		sys.exit(1)

	# run calls
	if run_calls(sys.argv[1], hist_fn, sys.argv[3]) != 0:
		sys.exit(1)
	sys.exit(0)
# end of main

