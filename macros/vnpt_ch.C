#ifndef VNPT
#define VNPT
#include <signal.h>
#include "Framework.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TH1I.h"
#include <ctime>
#include "../CMSSW_10_3_1_patch3/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;
using namespace std;
static const int nptbins = 12;
static const float ptbins[]={0.5, 0.6,  0.8,  1.0,  1.40,  1.8,
			      2.2,  2.8,  3.6,  4.6,  6.0,  7.0, 8.5};

float sub1[2] = {-1.,0};
float sub2[2] = {0,1.0};
int ROI[nptbins];
int ptloc[100];
TFile * fout = NULL;
TGraphErrors * vnptCalc(string flist, int mintrk, int maxtrk, int *Markers, int *Colors);

void vnpt_ch(){
  int Markers[2];
  int Colors[2];

  fout = new TFile("ch.root","recreate");
  Markers[0] = 20;
  Markers[1] = 21;
  Colors[0] = kRed;
  Colors[1] = kBlue;
  vnptCalc("filelists/PbPb_ch.dat", 15, 30, Markers, Colors);
}

TGraphErrors * vnptCalc(string flist, int minCent, int maxCent, int *Markers, int *Colors){
  std::clock_t start = 0;
  std::clock_t last = 0;
  int order = 2;
  Framework * frame = new Framework(flist);
  int iloc = 0;
  
 
  for(int j = 0; j<nptbins; j++) {
    ROI[j] = frame->SetROIRange(order, minCent,maxCent, sub1,sub2,  ptbins[j],ptbins[j+1]);
    frame->SetROIEP(ROI[j],HFp2,HFm2,trackmid2,HFm2,HFp2,trackmid2);
    if(ROI[j]<0) {cout<<"ROI ERROR"<<endl; return NULL;}
    ptloc[iloc++] = j;
  }
  cout<<"iloc: "<<iloc<<endl;
  int nfiles = 0;
  while(frame->AddFile()) {++nfiles;}
  cout<<"processed "<<nfiles<< " files"<<endl;
  double x[2][17];
  double xSub[2][17];
  double y[2][17];
  double ySub[2][17];
  double ey[2][17];
  double eySub[2][17];
  TGraphErrors * g[2];
  TGraphErrors * gSub[2];
  for(int isub = 0; isub<=1; isub++) {
    for(int i = 0; i<nptbins; i++) {
      x[isub][i] = frame->GetAvPt(i,0);
      xSub[isub][i] = x[isub][i];
      y[isub][i] = frame->GetVn(i,isub);
      ey[isub][i] = frame->GetVnErr(i,isub);
      ySub[isub][i] = frame->GetVnSubEvt(i,isub);
      eySub[isub][i] = frame->GetVnErrSubEvt(i,isub);
      cout<<isub<<"\t"<<i<<"\t"<<x[isub][i]<<"\t"<<y[isub][i]<<endl;
    }
    g[isub] = new TGraphErrors(nptbins,x[isub],y[isub],0,ey[isub]);
    gSub[isub] = new TGraphErrors(nptbins,xSub[isub],ySub[isub],0,eySub[isub]);
    g[isub]->SetName(Form("g_sub%d",isub+1));
    g[isub]->SetTitle(Form("g_sub%d",isub+1));
    g[isub]->SetMarkerStyle(Markers[isub]);
    g[isub]->SetMarkerColor(Colors[isub]);
    gSub[isub]->SetName(Form("gSub_sub%d",isub+1));
    gSub[isub]->SetTitle(Form("gSub_sub%d",isub+1));
    gSub[isub]->SetMarkerStyle(Markers[isub]);
    gSub[isub]->SetMarkerColor(Colors[isub]);
  }
  
  fout->cd();
  TDirectory * subdir;
  subdir = fout->mkdir(Form("PbPb_%d_%d",minCent,maxCent));
  subdir->cd();
  g[0]->Write("vn_NegEta");
  gSub[0]->Write("vn_NegEta_SubEvt");
  
  g[1]->Write("vn_PosEta");
  gSub[1]->Write("vn_PosEta_SubEvt");

  //spec->Write();
  
  return 0;
}

#endif
