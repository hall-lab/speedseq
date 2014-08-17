#ifndef __GENOTYPER_HH__
#define __GENOTYPER_HH__

#include "HisMaker.hh"

class Genotyper
{
private:
  HisMaker *_maker;
  TString  _file;
  int      _bin;
  TH1 *_hisSignal;
  TH1 *_hisDistr,*_hisDistrAll;
  TH1 *_hisDistr1000,*_hisDistr1000All;
  double _mean,   _sigma;
  double _meanAll,_sigmaAll;
  double _mean1000,   _sigma1000;
  double _mean1000All,_sigma1000All;

public:
  Genotyper(HisMaker *maker,TString file,int bin);
  void printGenotype(TString chrom,int start,int end,
		     bool useATcorr,bool useGCcorr);

private:
  double getReadCount(int start,int end);
};

#endif
