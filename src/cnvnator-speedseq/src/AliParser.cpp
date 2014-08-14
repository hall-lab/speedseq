#include "AliParser.hh"

AliParser::AliParser(string fileName,bool loadIndex) : sam(false),
						       bam(false),
						       stdin(false),
						       file(NULL),
						       index(NULL),
						       record(NULL),
						       fin(NULL),
						       samin(NULL),
						       cnames_(NULL),
						       clens_(NULL),
						       n_chr_(0),
						       flag_(0)
{
  int len = fileName.length();
  if (len == 0) {
    stdin = true;
  } else if (fileName.substr(len - 3,3) == "bam") {
    //cout<<"Assuming BAM file for "<<fileName<<endl;
    file = samopen(fileName.c_str(),"rb",NULL);
    if (!file) cerr<<"Can't open file '"<<fileName<<"'."<<endl;
    else {
      if (file->header->n_targets <= 0) {
	cerr<<"Header contains no sequences."<<endl;
	cerr<<"Aborting parsing."<<endl;
      } else {
	bam = true;
	record = new bam1_t();
      }
    }
    if (loadIndex) {
      index = bam_index_load(fileName.c_str());
      if (!index) cerr<<"Can't load index file for '"<<fileName<<"'."<<endl;
    }
  } else if (fileName.substr(len - 3,3) == "sam") {
    cout<<"Assuming SAM file ..."<<endl;
    file = samopen(fileName.c_str(),"r",NULL);
    if (!file) cerr<<"Can't open file '"<<fileName<<"'."<<endl;
    else {
      sam = true;
      if (file->header->n_targets <= 0) {
	cerr<<"Header contains no sequences in file '"<<fileName<<"'."<<endl;
	cerr<<"Parsing line by line ..."<<endl;
	samin = new ifstream(fileName.c_str());
	if (!samin->good()) {
	  cerr<<"Can't open file '"<<fileName<<"'."<<endl;
	  samin = NULL;
	}
      } else record = new bam1_t();
    }
  } else {
    fin = new ifstream(fileName.c_str());
    if (!fin->good()) {
      cerr<<"Can't open file '"<<fileName<<"'."<<endl;
      fin = NULL;
    }
  }
  if (file) { // sam or bam
      n_chr_ = file->header->n_targets;
      cnames_ = new string[n_chr_];
      clens_  = new int[n_chr_];
      for (int i = 0;i < n_chr_;i++) {
	clens_[i]  = file->header->target_len[i];
	cnames_[i] = file->header->target_name[i];
      }
  }
}

AliParser::~AliParser()
{
  if (fin)     delete fin;
  if (samin)   delete samin;
  if (file)    samclose(file);
  if (record)  delete record;
  if (index)   bam_index_destroy(index);
  if (cnames_) delete[] cnames_;
  if (clens_)  delete[] clens_;
}

int AliParser::parseRegion(bam_fetch_f callback, const string &chunk, void *data)
{
    int tid = 0, beg = 0, end = 0;
    if (bam_parse_region(this->file->header, chunk.c_str(), &tid, &beg, &end) < 0) {
        return -1;
    }
    bam_fetch(this->file->x.bam, this->index, tid, beg, end, data, callback);
    return 0;
}

bool AliParser::parseRecord()
{
  chr_index_ = -1;
  if (bam || (sam && !samin)) {
    if (samread(file,record) < 0) return false;
    flag_ = record->core.flag;
    bam1_core_t &core = record->core;
    chr_index_ = core.tid;
    start_     = core.pos + 1;
    end_       = bam_calend(&core,bam1_cigar(record));
    qual_      = core.qual;
    read_len_  = core.l_qseq;
    frg_len_   = core.isize;
    //if (frg_len_ < 0) frg_len_ = -frg_len_;
    if (chr_index_ >= 0 && chr_index_ < n_chr_)
      chr_ = cnames_[chr_index_];
    else chr_ = "?";
  } else if (stdin || samin) {
    if (samin) return parseSamLine(samin);
    else       return parseSamLine(&cin);
  } else {
    char c;
    chr_index_ = -1;
    if (!fin) return false;
    if (!fin->eof()) *fin>>chr_;      else return false;
    if (!fin->eof()) *fin>>start_;    else return false;
    if (!fin->eof()) *fin>>read_len_; else return false;
    if (!fin->eof()) *fin>>c;         else return false;
    end_ = start_ + read_len_ - 1;
    if (c == 'U' || c == 'S') qual_ = 100;
    else                      qual_ =   0;
    frg_len_ = 0;
  }
  return true;
}

bool AliParser::parseSamLine(istream *sin)
{
  const static int max = 4096;
  static char buf[max];
  (*sin).getline(buf,max);
  if ((*sin).eof()) return false;

  if (buf[0] == '@') { // Comment
    chr_index_ =  -1;
    chr_       =  "";
    start_     =   0;
    end_       =   0;
    read_len_  =   0;
    frg_len_   =   0;
    qual_      =   0;
    return true;
  }

  int i = 0;
  while (i < max && buf[i] != '\0' && buf[i] != '\t') i++; i++; // 1
  while (i < max && buf[i] != '\0' && buf[i] != '\t') i++; i++; // 2
  
  string tmp = "";
  while (i < max && buf[i] != '\0' && buf[i] != '\t') tmp += buf[i++]; i++;
  if (tmp == "*") chr_ = tmp;
  else {
    chr_ = "";
    int n = tmp.length();
    for (int i = 0;i < n;i++) {
      char c = tmp.at(i);
      if (!isalnum(c)) break;
      chr_ += c;
    }
  }

  tmp = "";
  while (i < max && buf[i] != '\0' && buf[i] != '\t') tmp += buf[i++]; i++;
  start_ = strtol(tmp.c_str(),NULL,10);
  
  tmp = "";
  while (i < max && buf[i] != '\0' && buf[i] != '\t') tmp += buf[i++]; i++;
  qual_ = strtol(tmp.c_str(),NULL,10);
  
  tmp = "";
  while (i < max && buf[i] != '\0' && buf[i] != '\t') tmp += buf[i++]; i++;
  //if (tmp == "*") return false;
  
  while (i < max && buf[i] != '\0' && buf[i] != '\t') i++; i++; // 7
  while (i < max && buf[i] != '\0' && buf[i] != '\t') i++; i++; // 8
  while (i < max && buf[i] != '\0' && buf[i] != '\t') i++; i++; // 9
  
  int len = 0;
  while (i < max && buf[i] != '\0' && buf[i] != '\t') { len++; i++; }
  end_ = start_ + len - 1;
  
  if (i >= max) return false;

  return true;
}

// callback for bam_fetch()  
static int fetch_func(const bam1_t *b, void *data)
{
  return 0;
}
int AliParser::scrollTo(string chr,int start)
{
  int chri = -1;
  for (int i = 0;i < n_chr_;i++)
    if (chr == cnames_[i]) {
      chri = i;
      break;
    }
  if (chri < 0) return -1;
  int read_len = 150;
  if (start < read_len) start = read_len;
  bam_fetch(file->x.bam,index,chri,start - read_len,start - read_len + 1,
	    NULL,fetch_func);
  while (parseRecord()) {
    if (chr_index_ < 0) continue;
    if (chr_index_ == chri) return chri;
    else                    return -1;
  }
  return -1;
}
