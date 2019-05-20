#ifndef CHECKFILES
#define CHECKFILES
#ifndef DEBUG
#define DEBUG
bool debug = false;
#endif
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include "../CMSSW_10_3_2/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#ifndef EFFFILE
#define EFFFILE
string efffile = "eff.root";
#endif
using namespace hi;
using namespace std;
  static const int nqxorder = 1;
  int qxorders[nqxorder]={2};

class CheckFiles{
public:
  CheckFiles(string filelist="filelist.dat");
  bool AddFile();
  int GetN(){return maxevents;}
  void GetEntry(int i){tr->GetEntry(i);}
  FILE * flist;
private:
  TRandom * ran;
  int maxevents;
  TFile * tf;
  TTree * tr;
  double centval;
  int ntrkval;
  int Noff;
  double vtx;
  double epang[NumEPNames];
  Double_t eporig[NumEPNames];
  Double_t qx[NumEPNames];
  Double_t qy[NumEPNames];
  Double_t q[NumEPNames];
  Double_t vn[NumEPNames];
  Double_t epmult[NumEPNames];
  Double_t sumw[NumEPNames];
  Double_t sumw2[NumEPNames];
  Double_t rescor[NumEPNames];
  Double_t rescorErr[NumEPNames];
  unsigned int runno_;
  TH2F * qxtrk[7]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
  TH2F * qytrk[7]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
  TH2F * qcnt=nullptr;
  TH2F * avpt=nullptr;
};
bool CheckFiles::AddFile(){
  char buf[120];
  if(fgets(buf,120,flist)==NULL) return false;
  buf[strlen(buf)-1]=0;
  string infile = buf;
  tf->Close();
  tf = new TFile(infile.data(),"read");
  tr = (TTree *) tf->Get("vnanalyzer/tree");
  maxevents = tr->GetEntries();
  cout<<maxevents<<"\tinfile:"<<infile<<":"<<endl;
  return true;
}
CheckFiles::CheckFiles(string filelist){
  ran = new TRandom();
  cout<<"open: "<<filelist<<endl;
  system("cat filelist.dat");
  flist = fopen(filelist.data(),"r");
  char buf[120];
  fgets(buf,120,flist);
  buf[strlen(buf)-1]=0;
  string infile = buf;
  tf = new TFile(infile.data(),"read");
  tr = (TTree *) tf->Get("vnanalyzer/tree");
  maxevents = tr->GetEntries();
  cout<<maxevents<<"\tinfile:"<<infile<<":"<<endl;
  while(AddFile());  
}
#endif
