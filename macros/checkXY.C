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
#include "TEllipse.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "Framework.h"
#include "../CMSSW_10_3_2/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;
using namespace std;
static const int maxFiles = 2000;
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
TH2D * hDiff[ncentbins][nptbins];
TH1D * diffx[ncentbins][nptbins];
TH1D * diffy[ncentbins][nptbins];
TH1D * runqdifx[ncentbins][nptbins];
TH1D * runqdify[ncentbins][nptbins];
TH1D * runqdifcnt[ncentbins][nptbins];
TH1I * hruns;
double runlist[500];
int nruns = 0;
int ROIp[ncentbins][nptbins];
int ROIm[ncentbins][nptbins];
int mapp[ncentbins][nptbins];
int mapm[ncentbins][nptbins];
int centloc[500];
int ptloc[500];
Framework * frame;

void CreatePlots(int icent, int ipt);

void checkXY(string inlist="filelist.dat"){
  char buf[120];
  TFile * frun = new TFile("runs.root","read");
  runs = (TH1I *) frun->Get("runs");
  runs->SetDirectory(0);
  for(int i = 1; i<=runs->GetNbinsX(); i++) {
    if(runs->GetBinContent(i)>0) {
      runlist[nruns++] = (int) runs->GetBinCenter(i);
      cout<<(int) runs->GetBinCenter(i)<<"\t"<<setw(12)<<setprecision(8)<<runs->GetBinContent(i)<<endl;
    }
  }
  runlist[nruns] = runlist[nruns-1]+1;
  hruns = new TH1I("hruns","hruns",nruns,runlist);

  double centval;
  int ntrkval;
  int Noff;
  double vtx;
  unsigned int runno_;
  TTree * tr;
  //  system(Form("cat %s",inlist.data()));
  frame = new Framework(inlist);
  frame->SetMinMult(2);
  frame->SetRuns(nruns,runlist);
  if(!frame->LoadOffsets("checkXY_offsets.root")) {
    cout<<"Failed to open offsets file"<<endl; }
  
  int iloc = 0;
  for(int i = 0; i<ncentbins; i++) {
    for(int j = 0; j<nptbins; j++) {
      ROIp[i][j]=frame->SetROIRange(Order,centbins[i],centbins[i+1],0.4, 0.8, ptbins[j], ptbins[j+1]);
      frame->SetROIEP(ROIp[i][j],HFm2,HFp2,trackmid2);
      mapp[i][j]=iloc;
      ptloc[iloc] = j;
      centloc[iloc++]=i;
      ROIm[i][j]=frame->SetROIRange(Order,centbins[i],centbins[i+1],-0.8, -0.4, ptbins[j], ptbins[j+1]);
      frame->SetROIEP(ROIm[i][j],HFp2,HFm2,trackmid2);
      mapm[i][j]=iloc;
      ptloc[iloc] = j;
      centloc[iloc++]=i;
      hDiff[i][j] = new TH2D(Form("diff_%d_%d_%03.1f_%03.1f",centbins[i],centbins[i+1],ptbins[j],ptbins[j+1]),
			     Form("diff_%d_%d_%03.1f_%03.1f",centbins[i],centbins[i+1],ptbins[j],ptbins[j+1]),200,-1,1,200,-1,1);
      hDiff[i][j]->SetDirectory(0);
      hDiff[i][j]->SetOption("colz");
      hDiff[i][j]->SetXTitle("(v_{2,x}^{obs,a}-v_{2,x}^{obs,b})/2");
      hDiff[i][j]->SetYTitle("(v_{2,y}^{obs,a}-v_{2,y}^{obs,b})/2");

      runqdifx[i][j] = new TH1D(Form("runqdifx_%d_%d_%03.1f_%03.1f",centbins[i],centbins[i+1],ptbins[j],ptbins[j+1]),
				Form("runqdifx_%d_%d_%03.1f_%03.1f",centbins[i],centbins[i+1],ptbins[j],ptbins[j+1]),nruns,runlist);
      runqdifx[i][j]->SetDirectory(0);
      runqdifx[i][j]->Sumw2();
      runqdify[i][j] = new TH1D(Form("runqdify_%d_%d_%03.1f_%03.1f",centbins[i],centbins[i+1],ptbins[j],ptbins[j+1]),
				Form("runqdify_%d_%d_%03.1f_%03.1f",centbins[i],centbins[i+1],ptbins[j],ptbins[j+1]),nruns,runlist);
      runqdify[i][j]->SetDirectory(0);
      runqdify[i][j]->Sumw2();
      runqdifcnt[i][j] = new TH1D(Form("runqdifcnt_%d_%d_%03.1f_%03.1f",centbins[i],centbins[i+1],ptbins[j],ptbins[j+1]),
				Form("runqdifcnt_%d_%d_%03.1f_%03.1f",centbins[i],centbins[i+1],ptbins[j],ptbins[j+1]),nruns,runlist);
      runqdifcnt[i][j]->SetDirectory(0);
      runqdifcnt[i][j]->Sumw2();
      


    }
  }
  frame->ShowAllROIRanges();
  int nfiles = 1;
  int count = 0;
  int partcount = 0;
 rep: for(int i = 0; i<frame->GetN(); i++) {
    signal(SIGINT,intHandler);
    if(keepRunning == 0) break;
    int runno = frame->GetRunno(i);
    int irun = hruns->FindBin(runno+0.1);
    frame->AddEvent(i);
    for(int k = 0; k<ncentbins; k++) {
      for(int j = 0; j<nptbins; j++) {
	if(frame->GetVnxEvt(mapp[k][j])<-1) continue;
	if(frame->GetVnxEvt(mapm[k][j])<-1) continue;
	if(frame->GetVnyEvt(mapp[k][j])<-1) continue;
	if(frame->GetVnyEvt(mapm[k][j])<-1) continue;
	double xdiffoff = frame->GetXdiff(mapp[k][j],runno);
	double ydiffoff = frame->GetYdiff(mapp[k][j],runno);
	double xdiff = (frame->GetVnxEvt(mapp[k][j])-frame->GetVnxEvt(mapm[k][j]))/2.;
	double ydiff = (frame->GetVnyEvt(mapp[k][j])-frame->GetVnyEvt(mapm[k][j]))/2.;
 	hDiff[k][j]->Fill(xdiff-xdiffoff ,ydiff-ydiffoff); 
	runqdifx[k][j]->Fill(runno,xdiff-xdiffoff);
	runqdify[k][j]->Fill(runno,ydiff-ydiffoff);
	runqdifcnt[k][j]->Fill(runno);
      }
    }
    ++count;
    ++partcount;
    if(partcount == 100000) {
      cout<<count<<endl;
      partcount = 0;
    }
  }
  if(nfiles++<maxFiles && keepRunning == 1 && frame->AddFile()) goto rep;
  TFile * fout = new TFile("checkXY.root","recreate");
  fout->cd();
  TDirectory * subdirs[ncentbins];
  TDirectory * subsubdirs[ncentbins][nptbins];
  for(int i = 0; i<ncentbins; i++) {
    subdirs[i] = fout->mkdir(Form("%d_%d",centbins[i],centbins[i+1]));
    subdirs[i]->cd();
    TH1D * sp = frame->GetSpectra(ROIp[0][i]);
    sp->Write();
    for(int j = 0; j<nptbins; j++) {
      subsubdirs[i][j] = subdirs[i]->mkdir(Form("%03.1f_%03.1f",ptbins[j],ptbins[j+1]));
    }
  }
   for(int i = 0; i<frame->GetNrange(); i++) {
     subsubdirs[centloc[i]][ptloc[i]]->cd();
     ((TH2D *)frame->Get2d(i))->Write();
     ((TH1D *)frame->GetMult(i))->Write();
   }


  for(int i = 0; i<ncentbins; i++) {
    for(int j = 0; j<nptbins; j++) {
      subsubdirs[i][j]->cd();
      hDiff[i][j]->Write("diff2d");
      runqdifx[i][j]->Write("qdifx");
      runqdify[i][j]->Write("qdify");
      runqdifcnt[i][j]->Write("qdifcnt");
      CreatePlots(i,j);
    }
  }
  for(int i = 0; i<frame->GetNrange(); i++) {
    frame->SaveFrameworkRuns(i,subsubdirs[centloc[i]][ptloc[i]]);
  }

}

void   CreatePlots(int icent, int ipt){

  string cname = Form("xy_%d_%d_%3.1f_%3.1f",centbins[icent],centbins[icent+1],ptbins[ipt],ptbins[ipt+1]);
  TCanvas * clow = new TCanvas(cname.data(),cname.data(),1300,1300);
  clow->Divide(2,2);
  clow->cd(1);

  frame->Get2d(mapm[icent][ipt])->Draw();
  gPad->SetGrid(1,1);
  TLine * l1a = new TLine(-1.4,0,1.4,0);
  l1a->SetLineWidth(2);
  l1a->Draw();
  TLine * l2a = new TLine(0,-1.4,0,1.4);
  l2a->SetLineWidth(2);
  l2a->Draw();
  TEllipse * el1a = new TEllipse(0,0,1,1);
  el1a->SetFillStyle(0);
  el1a->SetLineWidth(3);
  el1a->SetLineStyle(2);
  el1a->SetLineColor(kRed);
  el1a->Draw();
  clow->cd(2);
  frame->Get2d(mapp[icent][ipt])->Draw();
  gPad->SetGrid(1,1);
  TLine * l1b = new TLine(-1.4,0,1.4,0);
  l1b->SetLineWidth(2);
  l1b->Draw();
  TLine * l2b = new TLine(0,-1.4,0,1.4);
  l2b->SetLineWidth(2);
  l2b->Draw();
  TEllipse * el1b = new TEllipse(0,0,1,1);
  el1b->SetFillStyle(0);
  el1b->SetLineWidth(3);
  el1b->SetLineStyle(2);
  el1b->SetLineColor(kRed);
  el1b->Draw();
  clow->cd(3);
  hDiff[icent][ipt]->Draw();
  gPad->SetGrid(1,1);
  TLine * l1c = new TLine(-1.4,0,1.4,0);
  l1c->SetLineWidth(2);
  l1c->Draw();
  TLine * l2c = new TLine(0,-1.4,0,1.4);
  l2c->SetLineWidth(2);
  l2c->Draw();
  clow->cd(4);
  TPaveText * tp = new TPaveText(0.2,0.2,0.8,0.8);
  tp->SetBorderSize(0);
  tp->SetFillColor(kWhite);
  tp->AddText("2018 PbPb");
  tp->AddText(Form("%d - %d%c",centbins[icent],centbins[icent+1],'%'));
  tp->AddText(Form("%3.1f < p_{T} <  %3.1f GeV/c",ptbins[ipt],ptbins[ipt+1]));
  tp->Draw();
  clow->Print(Form("figures/checkXY/%s.pdf",clow->GetName()),"pdf");
  clow->Close();
  gSystem->ProcessEvents();
  delete clow;
  clow = 0;
  return;
}
