#ifndef __ALIPARSER__
#define __ALIPARSER__

// C/C++ includes
#include <iostream> 
#include <fstream> 
using namespace std; 

// Samtools includes
#include "sam.h"

class AliParser
{
private:
  bool         sam,bam,stdin;
  ifstream    *fin;
  ifstream    *samin;
  samfile_t   *file;
  bam_index_t *index;
  bam1_t      *record;
  string      *cnames_;
  int         *clens_;
  int          n_chr_;

public:
  int    numChrom() { return n_chr_; }
  string chromName(int i) { return (i >= 0 && i < n_chr_) ? cnames_[i] : ""; }
  int    chromLen(int i)  { return (i >= 0 && i < n_chr_) ? clens_[i] : 0; }

private: // Chromosome
  string chr_;
public:
  inline string getChromosome() { return chr_; }

private: // Chromosome index
  int chr_index_;
public:
  inline int getChromosomeIndex() { return chr_index_; }

private: // Mapping location start
  int start_;
public:
  inline int getStart() { return start_; }

private: // Mapping location end
  int end_;
public:
  inline int getEnd() { return end_; }

private: // Read length
  int read_len_;
public:
  inline int getReadLength() { return read_len_; }

private: // Fragment length
  int frg_len_;
public:
  inline int getFragmentLength() { return frg_len_; }

private: // Flag about read placement
  int qual_;
public:
  inline int getQuality() { return qual_; }
  inline bool isQ0() { return (qual_ <= 1); }

private: // Flag
  int flag_;
public:
  inline bool isUnmapped()     { return flag_ & 0x4; }
  inline bool isNextUnmapped() { return flag_ & 0x8; }
  inline bool isReversed()     { return flag_ & 0x10; }
  inline bool isNextReversed() { return flag_ & 0x20; }
  inline bool isDuplicate()    { return flag_ & 0x400; }

public:
  inline string getQueryName() { return (record) ? bam1_qname(record) : ""; }

public:
  AliParser(string fileName,bool loadIndex = false);
  ~AliParser();

  bool parseRecord();
  int  scrollTo(string chrom,int start);

private:
  bool parseSamLine(istream *sin);
    
public:
    int parseRegion(bam_fetch_f callback, const string &chunk, void *data);
};

#endif
