#include "TreeBuilder.hh"
#include "AliParser.hh"
#include <stdlib.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THashTable.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

typedef struct {
    string file;    // file name to parser
    string name;    // chrom name
    string cname;   // canonical chrom name
    uint32_t clen;  // chrom len
    uint64_t n_placed;  // number of alignments processed
    short *counts_p;    // counts
    short *counts_u;    // counts (unique)
    TH2 *his_frg_read;  // 2D histogram
} chrom_hist_data;


TreeBuilder::TreeBuilder(Genome *genome, string root_fn, int max_threads, bool forUnique)
{
    this->refGenome_ = genome;
    this->root_ = root_fn;
    this->max_threads_ = max_threads;
    this->forUnique_ = forUnique;
}


TreeBuilder::~TreeBuilder()
{
    for (int idx = 0; idx < counts_.size(); idx++) {
        chrom_hist_data *data = (chrom_hist_data *)counts_[idx];
        counts_[idx] = NULL;
        if (data->counts_u != NULL) {
            delete(data->counts_u);
        }
        delete(data->his_frg_read);
        delete(data->counts_p);
        delete(data);
    }
    counts_.clear();
}


int TreeBuilder::findIndex(string *arr, int n, string name)
{
    string can_name = Genome::makeCanonical(name);
    for (int i = 0; i < n; i++) {
        if (Genome::makeCanonical(arr[i]) == can_name) {
            return i;
        }
    }
    return -1;
}


int TreeBuilder::findIndex(const vector<string> &arr, string name)
{
    
    string can_name = Genome::makeCanonical(name);
    for (int i = 0; i < arr.size(); i++) {
        if (Genome::makeCanonical(arr[i]) == can_name) {
            return i;
        }
    }
    return -1;
}


void TreeBuilder::writeTreeForChromosome(string chrom, short *arr_p, short *arr_u, int len)
{
    // Creating a tree
    TFile file(root_.c_str(), "Update");
    if (file.IsZombie()) {
        cerr<<"Can't open/write to file '"<<root_<<"'."<<endl;
        return;
    }
    stringstream ss;
    ss<<chrom<<';'<<len;
    string description = ss.str();
    short rd_u,rd_p;
    int position;
    TTree *tree = new TTree(chrom.c_str(),description.c_str());
    tree->SetMaxTreeSize(20000000000); // ~20 Gb
    tree->Branch("position", &position, "position/I");
    tree->Branch("rd_unique",&rd_u,"rd_u/S");
    tree->Branch("rd_parity",&rd_p,"rd_p/S");
    // Filling the tree
    for (int i = 0;i < len;i++) {
        rd_u = (arr_u) ? arr_u[i] : 0;
        rd_p = (arr_p) ? arr_p[i] : 0;
        if (rd_u > 0 || rd_p > 0) {
            position = i;
            tree->Fill();
        }
    }
    
    // Writing the tree
    tree->Write(chrom.c_str(),TObject::kOverwrite);
    
    // Deleting the tree
    delete tree;
    file.Close();
}


int align_fetch(const bam1_t *align, void *data)
{
    AlignmentData x = AlignmentData(align);
    if (x.isUnmapped() || x.isDuplicate()) {
        return 0;
    }
    
    // update counts
    chrom_hist_data *hist_data = (chrom_hist_data *)data;
    int mid = abs(x.getStart() + x.getEnd()) >> 1;
    if (mid < 0 || mid > hist_data->clen) {
        return 0; // out of bounds
    }
    
    if (hist_data->counts_p[mid] + 1 > 0) {
        hist_data->counts_p[mid]++;
    }
    if (hist_data->counts_u && !x.isQ0()) {
        if (hist_data->counts_u[mid] + 1 > 0) {
            hist_data->counts_u[mid]++;
        }
    }
    hist_data->n_placed++;
    int frg_len = x.getFragmentLength();
    if (frg_len < 0) {
        hist_data->his_frg_read->Fill(x.getReadLength(),-frg_len);
    } else {
        hist_data->his_frg_read->Fill(x.getReadLength(), frg_len);
    }
    return 0;
}


#define THREAD_NOT_STARTED -1
#define THREAD_RUNNING 0
#define THREAD_DONE 1
#define THREAD_CLEAN 2

typedef struct {
    chrom_hist_data *data;
    char *state;
} thd_data;


void *tree_thread_run(void *data)
{
    thd_data *tdata = (thd_data *)data;
    AliParser parser = AliParser(tdata->data->file, true);
    // cout << "Processing BAM data corresponding to chrom " << tdata->data->cname << endl;
    parser.parseRegion(align_fetch, tdata->data->name, tdata->data);
    
    // clean up
    printf("%llu alignments were processed corresponding to chrom %s\n", tdata->data->n_placed, tdata->data->cname.c_str());
    *(tdata->state) = THREAD_DONE;
    delete(tdata);
    return NULL;
}


void TreeBuilder::buildParallel(const string &fn, const vector<pair<string, size_t> > &chroms)
{
    size_t ncs = chroms.size();
    vector<pthread_t> tid(ncs);
    vector<char> tstate(ncs);
    vector<chrom_hist_data *> tdata;
    
    // fill data for each chrom
    // note that each chrom will get a separate thread
    for (int x = 0; x < ncs; x++) {
        tstate[x] = THREAD_NOT_STARTED;
        
        chrom_hist_data *hist_data = new chrom_hist_data();
        hist_data->counts_p = (short *)calloc(chroms[x].second + 1, sizeof(short));
        if (forUnique_) {
            hist_data->counts_u = (short *)calloc(chroms[x].second + 1, sizeof(short));
        } else {
            hist_data->counts_u = NULL;
        }
        string hist_name = string("read_frg_len_") + chroms[x].first;
        hist_data->his_frg_read = new TH2I(hist_name.c_str(),"Read and fragment lengths", 300,0.5,300.5,3001,-0.5,3000.5);
        hist_data->file = string(fn);
        hist_data->name = chroms[x].first;
        hist_data->cname = Genome::makeCanonical(chroms[x].first);
        hist_data->clen = (uint32_t)chroms[x].second;
        hist_data->n_placed = 0;
        tdata.push_back(hist_data);
        counts_.push_back(hist_data);
        
        memset(&tid[x], 0x00, sizeof(pthread_t));
    }
    
    // run up to max_threads threads at once
    int curr_threads = 0;
    int thread_idx = 0;
    int finished = 0;
    while (true)
    {
        // create new threads if we have space
        while (curr_threads < max_threads_ && thread_idx < ncs) {
            thd_data *data = new thd_data(); // thread will clean this up
            data->data = tdata[thread_idx];
            data->state = &tstate[thread_idx];
            *(data->state) = THREAD_RUNNING;
            pthread_create(&tid[thread_idx], NULL, tree_thread_run, data);
            thread_idx++;
        }
        
        // check for the threads that have exited
        for (int idx = 0; idx < ncs; idx++) {
            if (tstate[idx] == THREAD_DONE) {
                void *foo = NULL;
                // don't start any new threads until we're sure that the others have exited
                (void)pthread_join(tid[idx], &foo);
                tstate[idx] = THREAD_CLEAN;
                finished++;
            }
        }
        if (finished == ncs) {
            break; // we've processed everything
        }
        curr_threads = thread_idx - finished;
        sleep(1); // don't peg the CPU
    }
}


int TreeBuilder::build(string *user_chroms, int chroms_len, string user_file)
{
    if (user_chroms == NULL) {
       chroms_len = 0;
    }
    
    // build a list of chromosome names and lengths
    // these are the chromosomes that we need to process
    AliParser *parser = new AliParser(user_file.c_str(), true);
    vector<pair<string, size_t> > chroms;
    if (parser->numChrom() == 0) {
        cout << "No chromosome / contig description given." << endl;
        return -1;
    }
    // use the parser's chroms as a basis
    // proc only the user-specified chroms list, or proc all chroms from the parser if the user didn't specify anything
    for (int c = 0; c < parser->numChrom(); c++) {
        string name = parser->chromName(c);
        if (chroms_len > 0 && findIndex(user_chroms, chroms_len, name) < 0) {
            continue;
        }
        chroms.push_back(pair<string, size_t>(name, parser->chromLen(c)));
    }
    this->buildParallel(user_file, chroms);
    return 0;
}


bool TreeBuilder::syncData()
{
    uint64_t n_placed = 0;
    TObjArray list = TObjArray(); // used to hold all of the individual histograms that we've accumulated
    for (int idx = 0; idx < counts_.size(); idx++) {
        chrom_hist_data *data = (chrom_hist_data *)counts_[idx];
        cout << "Filling and saving tree for '" << data->cname << "' ..." << endl;
        short *arru = NULL;
        if (data->counts_u != NULL) {
            arru = &data->counts_u[1];
        }
        short *arrp = &data->counts_p[1];
        writeTreeForChromosome(data->cname, arrp, arru, data->clen);
        list.Add(data->his_frg_read);
        n_placed += data->n_placed;
    }
    cout << "A total of " << n_placed << " aligned reads have been processed." << endl;
    
    // build histograms
    cout << "Writing histograms ... " << endl;
    static const int WIN = 2000;
    int win = WIN / 2;
    // we have to heap alloc these lest we get pages and pages of runtime warnings
    TH2 *his_at_aggr = new TH2I("his_at_aggr","AT aggregation", 51, 9.5, 60.5, 2*WIN + 1, -WIN - 0.5, WIN + 0.5);
    TH2 *his_frg_read = new TH2I("read_frg_len","Read and fragment lengths", 300, 0.5, 300.5, 3001, -0.5, 3000.5);
    TH3 *his_pair_pos = new TH3I("pair_pos","Pair position relative to AT run", 51, 9.5, 60.5,
                             2*win + 1, -win - 0.5, win + 0.5,
                             2*win + 1, -win - 0.5, win + 0.5);
    // merge all individual histograms into one
    his_frg_read->Merge(&list);
    
    TFile file(root_.c_str(), "Update");
    if (!file.IsZombie()) {
        // write histograms to file
        his_frg_read->Write(his_frg_read->GetName(), TObject::kOverwrite);
        his_at_aggr->Write(his_at_aggr->GetName(), TObject::kOverwrite);
        his_pair_pos->Write(his_pair_pos->GetName(), TObject::kOverwrite);
        file.Close();
    } else {
        cerr << "Can't open file '" << root_ << "'." << endl;
    }
    
    // clean up
    delete (his_at_aggr);
    delete (his_frg_read);
    delete (his_pair_pos);
    list.Clear();
    return true;
}