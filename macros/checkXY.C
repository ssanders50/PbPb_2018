#include <signal.h>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TH1I.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <iomanip>
#include "stdio.h"
#include "Framework.h"
#include "../CMSSW_10_3_2/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;
using namespace std;
static const int maxFiles = 4;
static volatile int keepRunning = 1;
void intHandler(int dummy){
  keepRunning = 0;
}
TFile * tf=NULL;
static const int nptbins = 17;
static const float ptbins[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0,
			      2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  7.0, 8.0, 10.};
static const int ncentbins = 13;
static const int centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,100};
static const int Order = 2;
FILE * flist;
TH1I * runs;
TH2D * hDiff[130];
void checkXY(string inlist="filelist.dat"){
  char buf[120];
  TFile * frun = new TFile("runs.root","read");
  runs = (TH1I *) frun->Get("runs");
  runs->SetDirectory(0);
  for(int i = 1; i<=runs->GetNbinsX(); i++) {
    if(runs->GetBinContent(i)>0) cout<<(int) runs->GetBinCenter(i)<<"\t"<<setw(12)<<setprecision(8)<<runs->GetBinContent(i)<<endl;
  }
  double centval;
  int ntrkval;
  int Noff;
  double vtx;
  unsigned int runno_;
  TTree * tr;
  int ROIp[nptbins][ncentbins];
  int ROIm[nptbins][ncentbins];
  //  system(Form("cat %s",inlist.data()));
  Framework * frame = new Framework(inlist);
  for(int i = 0; i<nptbins; i++) {
    for(int j = 0; j<ncentbins; j++) {
      ROIp[i][j]=frame->SetROIRange(Order,centbins[j],centbins[j+1],0., 0.8, ptbins[i], ptbins[i+1]);
      frame->SetROIEP(ROIp[i][j],HFm2,HFp2,trackmid2);
      ROIm[i][j]=frame->SetROIRange(Order,centbins[j],centbins[j+1],-0.8, 0., ptbins[i], ptbins[i+1]);
      frame->SetROIEP(ROIm[i][j],HFp2,HFm2,trackmid2);

    }
  }
  //frame->ShowAllROIRanges();
  int nfiles = 1;
  int count = 0;
  int partcount = 0;
 rep: for(int i = 0; i<frame->GetN(); i++) {
    signal(SIGINT,intHandler);
    if(keepRunning == 0) break;
    frame->AddEvent(i);
    ++count;
    ++partcount;
    if(partcount == 100000) {
      cout<<count<<endl;
      partcount = 0;
    }
  }
  if(nfiles++<maxFiles && keepRunning == 1 && frame->AddFile()) goto rep;
}
