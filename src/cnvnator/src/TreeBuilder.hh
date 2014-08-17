#ifndef __TreeBuilder__
#define __TreeBuilder__

#include <iostream>
#include <string>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THashTable.h>
#include <sam.h>
#include <sys/types.h>
#include <stdint.h>
#include "Genome.hh"
using std::string;

class AlignmentData
{
    public:
    AlignmentData(const bam1_t *data) { this->data = data; }
    
    inline string getQueryName() { return (data) ? bam1_qname(data) : ""; }
    inline bool isUnmapped()     { return data->core.flag & 0x4; }
    inline bool isNextUnmapped() { return data->core.flag & 0x8; }
    inline bool isReversed()     { return data->core.flag & 0x10; }
    inline bool isNextReversed() { return data->core.flag & 0x20; }
    inline bool isDuplicate()    { return data->core.flag & 0x400; }
    inline int32_t getChromosomeIndex() { return data->core.tid; }
    inline int32_t getStart() { return data->core.pos + 1; }
    inline int32_t getEnd() { return bam_calend(&data->core, bam1_cigar(data)); }
    inline int32_t getReadLength() { return data->core.l_qseq; }
    inline int32_t getFragmentLength() { return data->core.isize; }
    inline int32_t getQuality() { return data->core.qual; }
    inline bool isQ0() { return (data->core.qual <= 1); }
    
    private:
    const bam1_t *data;
};

class TreeBuilder
{
    public:
    TreeBuilder(Genome *genome, string root_fn, int max_threads = 32, bool forUnique = true);
    ~TreeBuilder();
    int build(string *user_chroms, int n_chroms, string user_files);
    bool syncData();
    
    private:
    void writeTreeForChromosome(string chrom, short *arr_p, short *arr_u, int len);
    int findIndex(string *arr, int n, string name);
    int findIndex(const vector<string> &arr, string name);
    void buildSerial(const string &fn, const vector<pair<string, size_t> > &chroms);
    void buildParallel(const string &fn, const vector<pair<string, size_t> > &chroms);
    
    private:
    Genome *refGenome_;
    string root_;
    int max_threads_;
    bool forUnique_;
    
    // histograms
    vector<void *> counts_;
};

#endif
