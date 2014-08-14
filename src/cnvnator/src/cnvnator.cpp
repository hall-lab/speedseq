// C/C++ includes
#include <iostream> 
#include <fstream> 
using namespace std; 

// Application includes
#include "AliParser.hh"
#include "HisMaker.hh"
#include "TreeBuilder.hh"
#define DEFAULT_MAX_THREADS 24
#define MAX_MAX_THREADS 48

int main(int argc,char *argv[])
{
  string usage = "\nCNVnator ";
#ifdef CNVNATOR_VERSION
  usage += CNVNATOR_VERSION;
#else
  usage += "v???";
#endif
  usage += "\n\nUsage:\n";
  usage += argv[0];
  usage += " -root out.root  [-genome name] [-chrom 1 2 ...] -tree  file1.bam ...\n";
  usage += argv[0];
  usage += " -root out.root  [-genome name] [-chrom 1 2 ...] [-nthreads threads] -ptree  file1.bam ...\n";
  usage += argv[0];
  usage += " -root out.root  [-genome name] [-chrom 1 2 ...] -merge file1.root ...\n";
  usage += argv[0];
  usage += " -root file.root [-genome name] [-chrom 1 2 ...] [-d dir] -his bin_size\n";
  usage += argv[0];
  usage += " -root file.root [-chrom 1 2 ...] -stat      bin_size\n";
  usage += argv[0];
  usage += " -root file.root                  -eval      bin_size\n";
  usage += argv[0];
  usage += " -root file.root [-chrom 1 2 ...] -partition bin_size [-ngc]\n";
  // usage += argv[0];
  //usage += " -root file.root [-chrom 1 2 ...] -spartition bin_size [-gc]\n";
  usage += argv[0];
  usage += " -root file.root [-chrom 1 2 ...] -call      bin_size [-ngc]\n";
  usage += argv[0];
  usage += " -root file.root -genotype bin_size [-ngc]\n";
  usage += argv[0];
  usage += " -root file.root -view     bin_size [-ngc]\n";
  usage += argv[0];
  usage += " -pe   file1.bam ... -qual val(20) -over val(0.8) [-f file]\n";
  usage += "\n";
  usage += "Valid genomes (-genome option) are: NCBI36, hg18, GRCh37, hg19\n";

  if (argc < 2) {
    cerr<<"Not enough parameters."<<endl;
    cerr<<usage<<endl;
    return 0;
  }

  static const int OPT_TREE       = 0x001;
  static const int OPT_MERGE      = 0x002;
  static const int OPT_HIS        = 0x004;
  static const int OPT_HISMERGE   = 0x008;
  static const int OPT_STAT       = 0x010;
  static const int OPT_PARTITION  = 0x020;
  static const int OPT_CALL       = 0x040;
  static const int OPT_VIEW       = 0x080;
  static const int OPT_GENOTYPE   = 0x100;
  static const int OPT_EVAL       = 0x200;
  static const int OPT_PE         = 0x400;
  static const int OPT_PTREE      = 0x800; // parallelized tree
  static const int OPT_NTHREADS   = 0x1000; // # of threads; useful with ptree

  static const int OPT_SPARTITION = 0x1000;
  static const int OPT_HIS_NEW    = 0x2000;
  static const int OPT_AGGREGATE  = 0x4000;
    

  // tree, merge, his, stat, partition, spartition, call, view, genotype
    int max_opts = 10000, n_opts = 0, max_threads = DEFAULT_MAX_THREADS;
    int opts[max_opts];
    int bins[max_opts];
    int gbin = 0;
  for (int i = 0;i < n_opts;i++) bins[i] = 0;
  bool useGCcorr = true,useATcorr = false;
  bool forUnique = false,relaxCalling = false;
  string out_root_file(""),call_file("");
  string chroms[1000],data_files[100000],root_files[100000] = {""},dir = ".";
  int n_chroms = 0,n_files = 0,n_root_files = 0,range = 128, qual = 20;
  double over = 0.8;
  Genome *genome = NULL;

  int index = 1;
  while (index < argc) {
    string option = argv[index++];
    if (option == "-tree"  || option == "-merge" || option == "-pe" || option == "-ptree") {
      if (option == "-tree")  opts[n_opts++] = OPT_TREE;
      if (option == "-merge") opts[n_opts++] = OPT_MERGE;
      if (option == "-pe")    opts[n_opts++] = OPT_PE;
      if (option == "-ptree") opts[n_opts++] = OPT_PTREE;
      while (index < argc && argv[index][0] != '-')
	if (strlen(argv[index++]) > 0) data_files[n_files++] = string(argv[index - 1]);
    } else if (option == "-his"       || option == "-his_new"    ||
	       option == "-hismerge"  ||
	       option == "-stat"      || option == "-eval"       ||
	       option == "-partition" || option == "-spartition" ||
	       option == "-call"      || option == "-view"       ||
	       option == "-genotype"  || option == "-aggregate") {
      int bs = 0;
      if (index < argc && argv[index][0] != '-') {
	TString tmp = argv[index++];
	if (!tmp.IsDigit()) {
	  cerr<<"Bin size must be integer for option '"<<option<<"'."<<endl;
	  cerr<<usage<<endl;
	  return 0;
	}
	bs = tmp.Atoi();
      }
      if (option == "-his")        opts[n_opts] = OPT_HIS;
      if (option == "-hismerge")   opts[n_opts] = OPT_HISMERGE;
      if (option == "-stat")       opts[n_opts] = OPT_STAT;
      if (option == "-partition")  opts[n_opts] = OPT_PARTITION;
      if (option == "-spartition") opts[n_opts] = OPT_SPARTITION;
      if (option == "-call")       opts[n_opts] = OPT_CALL;
      if (option == "-view")       opts[n_opts] = OPT_VIEW;
      if (option == "-genotype")   opts[n_opts] = OPT_GENOTYPE;
      if (option == "-his_new")    opts[n_opts] = OPT_HIS_NEW;
      if (option == "-eval")       opts[n_opts] = OPT_EVAL;
      if (option == "-aggregate")  opts[n_opts] = OPT_AGGREGATE;
      bins[n_opts++] = bs;
    } else if (option == "-root") {
      while (index < argc && argv[index][0] != '-')
	if (strlen(argv[index++]) > 0)
	  root_files[n_root_files++] = argv[index - 1];
      if (n_root_files == 0) {
	cerr<<"Please provide root-file name."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
    } else if (option == "-outroot") {
      if (index < argc && argv[index][0] != '-') out_root_file = argv[index++];
      else {
	cout<<"Please provide new root-file name."<<endl;
	return 0;
      }
    } else if (option == "-chrom") {
      while (index < argc && argv[index][0] != '-')
	chroms[n_chroms++] = argv[index++];
      if (n_chroms == 0) {
	cerr<<"Provide chromosome names."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
    } else if (option == "-ngc") {
      useGCcorr = false;
    } else if (option == "-at") {
      useATcorr = true;
    } else if (option == "-genome") {
      if (index < argc)	genome = Genome::get(argv[index++]);
    } else if (option == "-d") {
      if (index < argc && argv[index][0] != '-')
	dir = argv[index++];
      else cerr<<"No directory is given."<<endl;
    } else if (option == "-qual") {
      if (index >= argc || argv[index][0] == '-') {
	cerr<<"No quality value is provided."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      TString tmp = argv[index++];
      if (!tmp.IsDigit()) {
	cerr<<"Quality value must be integer."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      qual = tmp.Atoi();
    } else if (option == "-over") {
      if (index >= argc || argv[index][0] == '-') {
	cerr<<"No fraction of overlap is provided."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      TString tmp = argv[index++];
      if (!tmp.IsFloat()) {
	cerr<<"Fraction of overlap must be number."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      over = tmp.Atof();
    } else if (option == "-f") {
      if (index >= argc || argv[index][0] == '-') {
	cerr<<"No file name is provided."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      call_file = argv[index++];
    } else if (option == "-unique") {
      forUnique = true;
    } else if (option == "-range") {
      range = atoi(argv[index++]);
    } else if (option == "-relax") {
      relaxCalling = true;
    } else if (option == "-nthreads") {
        max_threads = atoi(argv[index++]);
        if (max_threads <= 0 || max_threads > MAX_MAX_THREADS) {
            cerr<<"Number of threads must be between 1 and "<<MAX_MAX_THREADS<<"."<<endl;
            cerr<<usage<<endl;
            return 0;
        }
    } else if (option[0] == '-') {
      cerr<<"Unknown option '"<<option<<"'.\n"<<endl;
    }
  }

  if (out_root_file.length() <= 0) out_root_file = root_files[0];
  if (out_root_file.length() <= 0)
    cerr<<"WARNING: no name of root-file provided."<<endl;

  for (int o = 0;o < n_opts;o++) {
    int option = opts[o];
    int bin = bins[o]; if (bin <= 0) bin = gbin;
    if (option == OPT_PTREE) { // tree
        TreeBuilder builder(genome, out_root_file, max_threads, forUnique);
        for (int idx = 0; idx < n_files; idx++) {
            builder.build(chroms, n_chroms, data_files[idx]);
        }
        (void)builder.syncData();
    }
    if (option == OPT_TREE) { // tree
        HisMaker maker(out_root_file,genome);
        maker.setDataDir(dir);
        maker.produceTrees(chroms,n_chroms,data_files,n_files,forUnique);
    }
    if (option == OPT_MERGE) { // merge
      HisMaker maker(out_root_file,genome);
      maker.mergeTrees(chroms,n_chroms,data_files,n_files);
    }
    if (option == OPT_HIS ||
	option == OPT_HISMERGE) { // his
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.setDataDir(dir);
      maker.produceHistograms(chroms,n_chroms,root_files,n_root_files,false);
      if (option == OPT_HISMERGE)
	maker.produceHistograms(chroms,n_chroms,root_files,n_root_files,true);
    }
    if (option == OPT_STAT) { // stat
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.stat(chroms,n_chroms,useATcorr);
    }
    if (option == OPT_PARTITION) { // partition
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.partition(chroms,n_chroms,false,useATcorr,useGCcorr,range);
    }
    if (option == OPT_CALL) { // call
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.callSVs(chroms,n_chroms,useATcorr,useGCcorr,relaxCalling);
    }
    if (option == OPT_VIEW) { // view
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      TApplication theApp("App",0,0);
      maker.view(root_files,n_root_files,useATcorr,useGCcorr);
      theApp.Run();
    }
    if (option == OPT_GENOTYPE) { // genotype
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      TApplication theApp("App",0,0);
      maker.genotype(root_files,n_root_files,useATcorr,useGCcorr);
      theApp.Run();
    }
    if (option == OPT_EVAL) { // eval
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.eval(root_files,n_root_files,useATcorr,useGCcorr);
    }
    if (option == OPT_PE) { // pe
      HisMaker maker("null",genome);
      if (call_file.length() > 0) 
	maker.pe_for_file(call_file,data_files,n_files,over,qual);
      else {
	TApplication theApp("App",0,0);
	maker.pe(data_files,n_files,over,qual);
	theApp.Run();
      }
    }
    if (option == OPT_SPARTITION) { // spartition
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.partition(chroms,n_chroms,true,useATcorr,useGCcorr,range);
    }
    if (option == OPT_HIS_NEW) { // his_new
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.setDataDir(dir);
      maker.produceHistogramsNew(chroms,n_chroms);
    }
    if (option == OPT_AGGREGATE) { // aggregate
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.setDataDir(dir);
      maker.aggregate(root_files,n_root_files,chroms,n_chroms);
    }
  }

  return 0;
}
