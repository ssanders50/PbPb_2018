#include <signal.h>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TH1I.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "stdio.h"
#include "../CMSSW_10_3_2/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;
using namespace std;
static volatile int keepRunning = 1;
void intHandler(int dummy){
  keepRunning = 0;
}
TFile * tf=NULL;
static int minRun = 326381;
static int maxRun = 327565;
int runcheck[400000];

FILE * flist;
TH1I * runs;
void scanFiles(string inlist="filelist.dat"){
  for(int i = 0; i<400000; i++) runcheck[i]=0;
  char buf[120];
  double centval;
  int ntrkval;
  int Noff;
  double vtx;
  unsigned int runno_;
  TTree * tr;
  runs = new TH1I("runs","runs",maxRun-minRun+1,minRun,maxRun);
  runs->SetDirectory(0);
  system(Form("cat %s",inlist.data()));
  flist = fopen(inlist.data(),"r");
  while(fgets(buf,120,flist)!=NULL) {
    buf[strlen(buf)-1]=0;
    cout<<buf<<endl;
    string infile = buf;
    if(tf!=NULL) tf->Close();
    tf = new TFile(infile.data(),"read");
    if(tf->IsZombie() || tf->TestBit(TFile::kRecovered) )                 {
      if(tf->IsZombie()) cout<<"ZOMBIE (removed):    " <<infile.data()<<endl; 
      if(tf->TestBit(TFile::kRecovered)) cout<<"Recovered (removed):    " <<infile.data()<<endl; 
      string remove = "rm "+infile;
      system(remove.data());
      continue;
    }
    tf->ResetErrno();
    tr = (TTree *) tf->Get("vnanalyzer/tree");
    int maxevents = tr->GetEntries();
    cout<<maxevents<<"\t"<<infile<<":"<<endl;
    tr->SetBranchAddress("Cent",&centval);
    tr->SetBranchAddress("NtrkOff",&Noff);
    tr->SetBranchAddress("ntrkflat",&ntrkval);
    tr->SetBranchAddress("Vtx",&vtx);
    tr->SetBranchAddress("Run",     &runno_);
    for(int i = 0; i<maxevents; i++) {
      tr->GetEntry(i);
      runs->Fill(runno_+0.5);
      ++runcheck[runno_];
    }
  }
  for(int i = 1; i<=runs->GetNbinsX(); i++) { 
    if(runs->GetBinContent(i)>0) cout<<runs->GetBinLowEdge(i)<<"\t"<<runs->GetBinContent(i)<<endl;
  }
  TFile * frun = new TFile("runs.root","recreate");
  frun->cd();
  runs->Write();
  frun->Close();
  for(int i = 0; i<400000; i++) {
    if(runcheck[i]!=0) cout<<"runcheck: "<<i<<"\t"<<runcheck[i]<<endl;
  }
}
