// Application includes
#include "Genome.hh"

Genome Genome::genomes[NGS] = {Genome("NCBI36"),Genome("GRCh37")};

Genome::Genome(string name) : n_chr_(0)
{
  string org_name = name;
  for(int i = 0;i < name.length();i++) name[i] = tolower(name[i]);
  if (name == "hg18" || name == "ncbi36") {
    gname_       = "NCBI36";
    other_gname_ = "hg18";
    n_chr_      = 24;
    cnames_[0]  =  "1"; clens_[0]  = 247249719;
    cnames_[1]  =  "2"; clens_[1]  = 242951149;
    cnames_[2]  =  "3"; clens_[2]  = 199501827;
    cnames_[3]  =  "4"; clens_[3]  = 191273063;
    cnames_[4]  =  "5"; clens_[4]  = 180857866;
    cnames_[5]  =  "6"; clens_[5]  = 170899992;
    cnames_[6]  =  "7"; clens_[6]  = 158821424;
    cnames_[7]  =  "8"; clens_[7]  = 146274826;
    cnames_[8]  =  "9"; clens_[8]  = 140273252;
    cnames_[9]  = "10"; clens_[9]  = 135374737;
    cnames_[10] = "11"; clens_[10] = 134452384;
    cnames_[11] = "12"; clens_[11] = 132349534;
    cnames_[12] = "13"; clens_[12] = 114142980;
    cnames_[13] = "14"; clens_[13] = 106368585;
    cnames_[14] = "15"; clens_[14] = 100338915;
    cnames_[15] = "16"; clens_[15] =  88827254;
    cnames_[16] = "17"; clens_[16] =  78774742;
    cnames_[17] = "18"; clens_[17] =  76117153;
    cnames_[18] = "19"; clens_[18] =  63811651;
    cnames_[19] = "20"; clens_[19] =  62435964;
    cnames_[20] = "21"; clens_[20] =  46944323;
    cnames_[21] = "22"; clens_[21] =  49691432;
    cnames_[22] =  "X"; clens_[22] = 154913754;
    cnames_[23] =  "Y"; clens_[23] =  57772954;
  } else if (name == "hg19" || name == "grch37") {
    gname_       = "GRCh37";
    other_gname_ = "hg19";
    n_chr_      = 24;
    cnames_[0]  =  "1"; clens_[0]  = 249250621;
    cnames_[1]  =  "2"; clens_[1]  = 243199373;
    cnames_[2]  =  "3"; clens_[2]  = 198022430;
    cnames_[3]  =  "4"; clens_[3]  = 191154276;
    cnames_[4]  =  "5"; clens_[4]  = 180915260;
    cnames_[5]  =  "6"; clens_[5]  = 171115067;
    cnames_[6]  =  "7"; clens_[6]  = 159138663;
    cnames_[7]  =  "8"; clens_[7]  = 146364022;
    cnames_[8]  =  "9"; clens_[8]  = 141213431;
    cnames_[9]  = "10"; clens_[9]  = 135534747;
    cnames_[10] = "11"; clens_[10] = 135006516;
    cnames_[11] = "12"; clens_[11] = 133851895;
    cnames_[12] = "13"; clens_[12] = 115169878;
    cnames_[13] = "14"; clens_[13] = 107349540;
    cnames_[14] = "15"; clens_[14] = 102531392;
    cnames_[15] = "16"; clens_[15] =  90354753;
    cnames_[16] = "17"; clens_[16] =  81195210;
    cnames_[17] = "18"; clens_[17] =  78077248;
    cnames_[18] = "19"; clens_[18] =  59128983;
    cnames_[19] = "20"; clens_[19] =  63025520;
    cnames_[20] = "21"; clens_[20] =  48129895;
    cnames_[21] = "22"; clens_[21] =  51304566;
    cnames_[22] =  "X"; clens_[22] = 155270560;
    cnames_[23] =  "Y"; clens_[23] =  59373566;
  } else {
    cerr<<"Unknown genome '"<<org_name<<"'."<<endl;
  }
}

int Genome::getChromosomeIndex(string chr)
{
  for (int c = 0;c < n_chr_;c++)
    if (chr == cnames_[c]) return c;
  return -1;
}

Genome *Genome::get(string name)
{
  for(int i = 0;i < name.length();i++) name[i] = tolower(name[i]);
  int ind = 0;
  while (ind < NGS) {
    string tmp = genomes[ind].gname_;
    for(int i = 0;i < tmp.length();i++) tmp[i] = tolower(tmp[i]);
    if (tmp == name) return &genomes[ind];
    tmp = genomes[ind].other_gname_;
    for(int i = 0;i < tmp.length();i++) tmp[i] = tolower(tmp[i]);
    if (tmp == name) return &genomes[ind];
    ind++;
  }
  return NULL;
}

string Genome::makeCanonical(string name)
{
  string tmp = "";
  for (int i = 0;i < name.length();i++) tmp += toupper(name[i]);
  string ret = "";
  if (tmp.substr(0,3)  == "CHR")             ret += name.substr(3);
  else if (tmp.substr(0,5)  == "CHROM")      ret += name.substr(5);
  else if (tmp.substr(0,10) == "CHROMOSOME") ret += name.substr(10);
  else if (name.length() < 3)                ret += name;
  else                                       ret  = name;
  if (ret == "chrx") ret = "chrX";
  if (ret == "chry") ret = "chrY";

  return ret;
}

