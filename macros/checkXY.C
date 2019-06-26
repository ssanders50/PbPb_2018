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
static const int maxFiles = 5;
static volatile int keepRunning = 1;
void intHandler(int dummy){
  keepRunning = 0;
}
TFile * tf=NULL;
static const int nptbins = 15;
static const float ptbinsMin[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0, 2.5,  3.0,  3.5, 0.3, 0.6, 1.0};
static const float ptbinsMax[]={0.4, 0.5, 0.6,  0.8,  1.0, 1.25,  1.50,   2.0,  2.5, 3.0,  3.5,  4.0, 3.0, 3.0, 3.0};
static const int ncentbins = 13;
static const int centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,100};
static const int Order = 2;
FILE * flist;
TH1I * runs;
TH2D * hDiff[ncentbins][nptbins];
TH2D * hDiffRef[ncentbins][nptbins];
TH2D * hResp[ncentbins][nptbins];
TH2D * hvncorr[ncentbins][nptbins];
TH2D * hangcorr[ncentbins][nptbins];
TH2D * hangcorrHFpA[ncentbins][nptbins];
TH2D * hangcorrHFpB[ncentbins][nptbins];
TH1D * hvn[ncentbins][nptbins];
TH1D * diffX[ncentbins][nptbins];
TH1D * diffY[ncentbins][nptbins];
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
      ROIp[i][j]=frame->SetROIRange(Order,centbins[i],centbins[i+1],0.0, 1.2, ptbinsMin[j], ptbinsMax[j]);
      frame->SetROIEP(ROIp[i][j],HFm2,HFp2,trackmid2);
      mapp[i][j]=iloc;
      ptloc[iloc] = j;
      centloc[iloc++]=i;
      ROIm[i][j]=frame->SetROIRange(Order,centbins[i],centbins[i+1],-1.2, 0.0, ptbinsMin[j], ptbinsMax[j]);
      frame->SetROIEP(ROIm[i][j],HFp2,HFm2,trackmid2);
      mapm[i][j]=iloc;
      ptloc[iloc] = j;
      centloc[iloc++]=i;
      hDiff[i][j] = new TH2D(Form("diff_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("diff_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),200,-1,1,200,-1,1);
      hDiff[i][j]->SetDirectory(0);
      hDiff[i][j]->SetOption("colz");
      hDiff[i][j]->SetXTitle("(v_{2,x}^{obs,a}-v_{2,x}^{obs,b})/2");
      hDiff[i][j]->SetYTitle("(v_{2,y}^{obs,a}-v_{2,y}^{obs,b})/2");

      hResp[i][j] = new TH2D(Form("resp_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("resp_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),100,0,1.0,100,0,1.0);
      hResp[i][j]->SetDirectory(0);
      hResp[i][j]->SetOption("colz");
      hResp[i][j]->SetXTitle("v_{2}^{obs}");
      hResp[i][j]->SetYTitle("v_{2}^{obs}");

      hvn[i][j] = new TH1D(Form("hvn_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			   Form("hvn_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),100,0,1);
      hvn[i][j]->SetDirectory(0);
      hvn[i][j]->SetOption("colz");
      hvn[i][j]->SetXTitle("(v_{2}^{a}+v_{2}^{b})/2");
      hvn[i][j]->SetYTitle("Counts");

      hvncorr[i][j] = new TH2D(Form("hvncorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hvncorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),100,0,1,100,0,1);
      hvncorr[i][j]->SetDirectory(0);
      hvncorr[i][j]->SetOption("colz");
      hvncorr[i][j]->SetXTitle("v_{2}^{a +}");
      hvncorr[i][j]->SetYTitle("v_{2}^{b -}");

      hangcorr[i][j] = new TH2D(Form("hangcorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				100,-2,2,100,-2,2);
      hangcorr[i][j]->SetDirectory(0);
      hangcorr[i][j]->SetOption("colz");
      hangcorr[i][j]->SetXTitle("#Psi^{b -}");
      hangcorr[i][j]->SetYTitle("#Psi^{a +}");


      hangcorrHFpA[i][j] = new TH2D(Form("hangcorrHFpA_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpA_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				100,-2,2,100,-2,2);
      hangcorrHFpA[i][j]->SetDirectory(0);
      hangcorrHFpA[i][j]->SetOption("colz");
      hangcorrHFpA[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpA[i][j]->SetYTitle("#Psi^{a +}");

      hangcorrHFpB[i][j] = new TH2D(Form("hangcorrHFpB_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpB_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				100,-2,2,100,-2,2);
      hangcorrHFpB[i][j]->SetDirectory(0);
      hangcorrHFpB[i][j]->SetOption("colz");
      hangcorrHFpB[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpB[i][j]->SetYTitle("#Psi^{b -}");


      diffX[i][j] = new TH1D(Form("diffX_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffX_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),400,-1.0,1.0);
      diffX[i][j]->SetDirectory(0);
      diffX[i][j]->SetOption("colz");
      diffX[i][j]->SetXTitle("(v_{2,x}^{obs,a}-v_{2,x}^{obs,b})/2");
      diffX[i][j]->SetYTitle("Counts");

      diffY[i][j] = new TH1D(Form("diffY_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffY_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),400,-1.0,1.0);
      diffY[i][j]->SetDirectory(0);
      diffY[i][j]->SetOption("colz");
      diffY[i][j]->SetXTitle("(v_{2,y}^{obs,a}-v_{2,y}^{obs,b})/2");
      diffY[i][j]->SetYTitle("Counts");


      runqdifx[i][j] = new TH1D(Form("runqdifx_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdifx_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdifx[i][j]->SetDirectory(0);
      runqdifx[i][j]->Sumw2();
      runqdify[i][j] = new TH1D(Form("runqdify_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdify_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdify[i][j]->SetDirectory(0);
      runqdify[i][j]->Sumw2();
      runqdifcnt[i][j] = new TH1D(Form("runqdifcnt_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdifcnt_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
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

	double vnxp = frame->GetVnxEvt(mapp[k][j]);
	double vnxm = frame->GetVnxEvt(mapm[k][j]);
	double vnx = (vnxp+vnxm)/2.;

	double vnyp = frame->GetVnyEvt(mapp[k][j]);
	double vnym = frame->GetVnyEvt(mapm[k][j]);
	double vny = (vnyp+vnym)/2.;

	double vnp = sqrt(pow(vnxp,2)+pow(vnyp,2));
	double vnm = sqrt(pow(vnxm,2)+pow(vnym,2));
	double vn = sqrt(pow(vnx,2)+pow(vny,2));
	hvn[k][j]->Fill(vn);
	hvncorr[k][j]->Fill(vnm,vnp);
	hangcorr[k][j]->Fill(frame->GetAng(mapm[k][j]),frame->GetAng(mapp[k][j]));
	hangcorrHFpA[k][j]->Fill(frame->GetAngHFp(mapp[k][j]),frame->GetAng(mapp[k][j]));
	hangcorrHFpB[k][j]->Fill(frame->GetAngHFp(mapp[k][j]),frame->GetAng(mapm[k][j]));
	double xdiff = (vnxp-vnxm)/2. - xdiffoff;
	double ydiff = (vnyp-vnym)/2. - ydiffoff;
 	hDiff[k][j]->Fill(xdiff,ydiff); 

	double vnobsx = 0;
	double vnobsy = 0;
	double vnobs = 0;
	frame->Smear(mapp[k][j],vnx,vny,vnobsx,vnobsy,vnobs);
	

	diffX[k][j]->Fill(xdiff,vnobs);
	diffY[k][j]->Fill(ydiff,vnobs);

	hResp[k][j]->Fill(vnobs,vn);
	runqdifx[k][j]->Fill(runno,xdiff);
	runqdify[k][j]->Fill(runno,ydiff);
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
      subsubdirs[i][j] = subdirs[i]->mkdir(Form("%04.2f_%04.2f",ptbinsMin[j],ptbinsMax[j]));
    }
  }
   for(int i = 0; i<frame->GetNrange(); i++) {
     subsubdirs[centloc[i]][ptloc[i]]->cd();
     ((TH2D *)frame->Get2d(i))->Write();
     ((TH1D *)frame->Get1d(i))->Write();
     TH1D * tmp = (TH1D *) frame->Get1dMult(i);
     tmp->Divide((TH1D *)frame->Get1d(i));
     tmp->Write();
     ((TH1D *)frame->GetMult(i))->Write();

   }


  for(int i = 0; i<ncentbins; i++) {
    for(int j = 0; j<nptbins; j++) {
      subsubdirs[i][j]->cd();
      hvn[i][j]->Write("vn");
      hvncorr[i][j]->Write("vncorr");
      hangcorr[i][j]->Write("angcorr");
      hangcorrHFpA[i][j]->Write("angcorr_HFpA");
      hangcorrHFpB[i][j]->Write("angcorr_HFpB");
      hDiff[i][j]->Write("diff2d");
      hResp[i][j]->Write("resp2d");
      diffX[i][j]->Write("diffX");
      diffY[i][j]->Write("diffY");
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

  string cname = Form("xy_%d_%d_%3.1f_%3.1f",centbins[icent],centbins[icent+1],ptbinsMin[ipt],ptbinsMax[ipt]);
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
  tp->AddText(Form("%3.1f < p_{T} <  %3.1f GeV/c",ptbinsMin[ipt],ptbinsMax[ipt]));
  tp->Draw();
  clow->Print(Form("figures/checkXY/%s.pdf",clow->GetName()),"pdf");
  clow->Close();
  gSystem->ProcessEvents();
  delete clow;
  clow = 0;
  return;
}
