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
#include "../CMSSW_10_3_3_patch1/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;
using namespace std;
static const int maxFiles = 5;
static volatile int keepRunning = 1;
void intHandler(int dummy){
  keepRunning = 0;
}
TFile * tf=NULL;
static const int nptbins = 15;
static const int nhistbins = 100;
static const double maxvn = 0.6;
static const float ptbinsMin[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0, 2.5,  3.0,  3.5, 0.3, 0.6, 1.0};
static const float ptbinsMax[]={0.4, 0.5, 0.6,  0.8,  1.0, 1.25,  1.50,   2.0,  2.5, 3.0,  3.5,  4.0, 3.0, 3.0, 3.0};
static const int ncentbins = 13;
static const int centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,100};
static const int Order = 2;
FILE * flist;
TH1I * runs;
TH2D * hDiff[ncentbins][nptbins];
TH2D * hDiff_p[ncentbins][nptbins];
TH2D * hDiff_m[ncentbins][nptbins];
TH2D * hDiffRef[ncentbins][nptbins];
TH2D * hDiffRef_p[ncentbins][nptbins];
TH2D * hDiffRef_m[ncentbins][nptbins];
TH2D * hResp[ncentbins][nptbins];
TH2D * hResp_p[ncentbins][nptbins];
TH2D * hResp_m[ncentbins][nptbins];
TH2D * hvncorr[ncentbins][nptbins];
TH2D * hvncorr_p[ncentbins][nptbins];
TH2D * hvncorr_m[ncentbins][nptbins];
TH2D * hangcorr[ncentbins][nptbins];
TH2D * hangcorr_p[ncentbins][nptbins];
TH2D * hangcorr_m[ncentbins][nptbins];
TH2D * hangcorrHFpA[ncentbins][nptbins];
TH2D * hangcorrHFpA_p[ncentbins][nptbins];
TH2D * hangcorrHFpA_m[ncentbins][nptbins];
TH2D * hangcorrHFpB[ncentbins][nptbins];
TH2D * hangcorrHFpB_p[ncentbins][nptbins];
TH2D * hangcorrHFpB_m[ncentbins][nptbins];
TH1D * hvn[ncentbins][nptbins]={0};
TH1D * hvn_p[ncentbins][nptbins]={0};
TH1D * hvn_m[ncentbins][nptbins]={0};
TH1D * diffX[ncentbins][nptbins];
TH1D * diffX_p[ncentbins][nptbins];
TH1D * diffX_m[ncentbins][nptbins];
TH1D * diffY[ncentbins][nptbins];
TH1D * diffY_p[ncentbins][nptbins];
TH1D * diffY_m[ncentbins][nptbins];
TH1D * runqdifx[ncentbins][nptbins];
TH1D * runqdifx_p[ncentbins][nptbins];
TH1D * runqdifx_m[ncentbins][nptbins];
TH1D * runqdify[ncentbins][nptbins];
TH1D * runqdify_p[ncentbins][nptbins];
TH1D * runqdify_m[ncentbins][nptbins];
TH1D * runqdifcnt[ncentbins][nptbins];
TH1D * runqdifcnt_p[ncentbins][nptbins];
TH1D * runqdifcnt_m[ncentbins][nptbins];
TH1D * hvnpt_pos_eta[ncentbins];
TH1D * hvnpt_pos_eta_p[ncentbins];
TH1D * hvnpt_pos_eta_m[ncentbins];
TH1D * hvnpt_neg_eta[ncentbins];
TH1D * hvnpt_neg_eta_p[ncentbins];
TH1D * hvnpt_neg_eta_m[ncentbins];
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

void CreatePlots(int icent, int ipt, string suffix = "");

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
  frame = new Framework(inlist,true, nhistbins, maxvn);
  frame->SetMinMult(2);
  frame->SetRuns(nruns,runlist);
  if(!frame->LoadOffsets("checkXY_offsets.root")) {
    cout<<"Failed to open offsets file"<<endl; }
  if(!frame->LoadOffsets_p("checkXY_p_offsets.root")) {
    cout<<"Failed to open offsets_p file"<<endl; }
  if(!frame->LoadOffsets_m("checkXY_m_offsets.root")) {
    cout<<"Failed to open offsets_m file"<<endl; }
  
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
			     Form("diff_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,-1,1,nhistbins,-1,1);
      hDiff[i][j]->SetDirectory(0);
      hDiff[i][j]->Sumw2();
      hDiff[i][j]->SetOption("colz");
      hDiff[i][j]->SetXTitle(Form("(v_{%d,x}^{obs,a}-v_{%d,x}^{obs,b})/2",Order,Order));
      hDiff[i][j]->SetYTitle(Form("(v_{%d,y}^{obs,a}-v_{%d,y}^{obs,b})/2",Order,Order));
      hDiff_p[i][j] = new TH2D(Form("diff_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("diff_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,-1,1,nhistbins,-1,1);
      hDiff_p[i][j]->SetDirectory(0);
      hDiff_p[i][j]->Sumw2();
      hDiff_p[i][j]->SetOption("colz");
      hDiff_p[i][j]->SetXTitle(Form("(v_{%d,x}^{obs,a}-v_{%d,x}^{obs,b})/2",Order,Order));
      hDiff_p[i][j]->SetYTitle(Form("(v_{%d,y}^{obs,a}-v_{%d,y}^{obs,b})/2",Order,Order));
      hDiff_m[i][j] = new TH2D(Form("diff_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("diff_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,-1,1,nhistbins,-1,1);
      hDiff_m[i][j]->SetDirectory(0);
      hDiff_m[i][j]->Sumw2();
      hDiff_m[i][j]->SetOption("colz");
      hDiff_m[i][j]->SetXTitle(Form("(v_{%d,x}^{obs,a}-v_{%d,x}^{obs,b})/2",Order,Order));
      hDiff_m[i][j]->SetYTitle(Form("(v_{%d,y}^{obs,a}-v_{%d,y}^{obs,b})/2",Order,Order));

      hResp[i][j] = new TH2D(Form("resp_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("resp_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn,nhistbins,0,maxvn);
      hResp[i][j]->SetDirectory(0);
      hResp[i][j]->Sumw2();
      hResp[i][j]->SetOption("colz");
      hResp[i][j]->SetXTitle(Form("v_{%d}^{obs}",Order));
      hResp[i][j]->SetYTitle(Form("v_{%d}^{obs}",Order));
      hResp_p[i][j] = new TH2D(Form("resp_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("resp_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn,nhistbins,0,maxvn);
      hResp_p[i][j]->SetDirectory(0);
      hResp_p[i][j]->Sumw2();
      hResp_p[i][j]->SetOption("colz");
      hResp_p[i][j]->SetXTitle(Form("v_{%d}^{obs}",Order));
      hResp_p[i][j]->SetYTitle(Form("v_{%d}^{obs}",Order));
      hResp_m[i][j] = new TH2D(Form("resp_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("resp_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn,nhistbins,0,maxvn);
      hResp_m[i][j]->SetDirectory(0);
      hResp_m[i][j]->Sumw2();
      hResp_m[i][j]->SetOption("colz");
      hResp_m[i][j]->SetXTitle(Form("v_{%d}^{obs}",Order));
      hResp_m[i][j]->SetYTitle(Form("v_{%d}^{obs}",Order));

      // hvn[i][j] = new TH1D(Form("hvn_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
      // 			   Form("hvn_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn);
      // hvn[i][j]->SetDirectory(0);
      // hvn[i][j]->Sumw2();
      // hvn[i][j]->SetOption("colz");
      // hvn[i][j]->SetXTitle(Form("(v_{%d}^{a}+v_{%d}^{b})/2",Order,Order));
      // hvn[i][j]->SetYTitle("Counts");
      // hvn_p[i][j] = new TH1D(Form("hvn_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
      // 			   Form("hvn_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn);
      // hvn_p[i][j]->SetDirectory(0);
      // hvn_p[i][j]->Sumw2();
      // hvn_p[i][j]->SetOption("colz");
      // hvn_p[i][j]->SetXTitle(Form("(v_{%d}^{a}+v_{%d}^{b})/2",Order,Order));
      // hvn_p[i][j]->SetYTitle("Counts");
      // hvn_m[i][j] = new TH1D(Form("hvn_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
      // 			   Form("hvn_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn);
      // hvn_m[i][j]->SetDirectory(0);
      // hvn_m[i][j]->Sumw2();
      // hvn_m[i][j]->SetOption("colz");
      // hvn_m[i][j]->SetXTitle(Form("(v_{%d}^{a}+v_{%d}^{b})/2",Order,Order));
      // hvn_m[i][j]->SetYTitle("Counts");

      hvnpt_pos_eta[i] = new TH1D(Form("hvnpt_pos_eta_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_pos_eta_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_pos_eta[i]->SetDirectory(0);
      hvnpt_pos_eta[i]->Sumw2();
      hvnpt_pos_eta[i]->SetOption("colz");
      hvnpt_pos_eta[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_pos_eta[i]->SetYTitle(Form("v_{%d} (+#eta)",Order));

      hvnpt_pos_eta_p[i] = new TH1D(Form("hvnpt_pos_eta_p_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_pos_eta_p_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_pos_eta_p[i]->SetDirectory(0);
      hvnpt_pos_eta_p[i]->Sumw2();
      hvnpt_pos_eta_p[i]->SetOption("colz");
      hvnpt_pos_eta_p[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_pos_eta_p[i]->SetYTitle(Form("v_{%d} (+#eta)",Order));

      hvnpt_pos_eta_m[i] = new TH1D(Form("hvnpt_pos_eta_m_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_pos_eta_m_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_pos_eta_m[i]->SetDirectory(0);
      hvnpt_pos_eta_m[i]->Sumw2();
      hvnpt_pos_eta_m[i]->SetOption("colz");
      hvnpt_pos_eta_m[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_pos_eta_m[i]->SetYTitle(Form("v_{%d} (+#eta)",Order));

      hvnpt_neg_eta[i] = new TH1D(Form("hvnpt_neg_eta_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_neg_eta_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_neg_eta[i]->SetDirectory(0);
      hvnpt_neg_eta[i]->Sumw2();
      hvnpt_neg_eta[i]->SetOption("colz");
      hvnpt_neg_eta[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_neg_eta[i]->SetYTitle(Form("v_{%d} (-#eta)",Order));

      hvnpt_neg_eta_p[i] = new TH1D(Form("hvnpt_neg_eta_p_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_neg_eta_p_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_neg_eta_p[i]->SetDirectory(0);
      hvnpt_neg_eta_p[i]->Sumw2();
      hvnpt_neg_eta_p[i]->SetOption("colz");
      hvnpt_neg_eta_p[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_neg_eta_p[i]->SetYTitle(Form("v_{%d} (-#eta)",Order));

      hvnpt_neg_eta_m[i] = new TH1D(Form("hvnpt_neg_eta_m_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_neg_eta_m_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_neg_eta_m[i]->SetDirectory(0);
      hvnpt_neg_eta_m[i]->Sumw2();
      hvnpt_neg_eta_m[i]->SetOption("colz");
      hvnpt_neg_eta_m[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_neg_eta_m[i]->SetYTitle(Form("v_{%d} (-#eta)",Order));

      hvncorr[i][j] = new TH2D(Form("hvncorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hvncorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn,nhistbins,0,maxvn);
      hvncorr[i][j]->SetDirectory(0);
      hvncorr[i][j]->Sumw2();
      hvncorr[i][j]->SetOption("colz");
      hvncorr[i][j]->SetXTitle("v_{2}^{a +}");
      hvncorr[i][j]->SetYTitle("v_{2}^{b -}");
      hvncorr_p[i][j] = new TH2D(Form("hvncorr_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hvncorr_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn,nhistbins,0,maxvn);
      hvncorr_p[i][j]->SetDirectory(0);
      hvncorr_p[i][j]->Sumw2();
      hvncorr_p[i][j]->SetOption("colz");
      hvncorr_p[i][j]->SetXTitle("v_{2}^{a +}");
      hvncorr_p[i][j]->SetYTitle("v_{2}^{b -}");
      hvncorr_m[i][j] = new TH2D(Form("hvncorr_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hvncorr_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn,nhistbins,0,maxvn);
      hvncorr_m[i][j]->SetDirectory(0);
      hvncorr_m[i][j]->Sumw2();
      hvncorr_m[i][j]->SetOption("colz");
      hvncorr_m[i][j]->SetXTitle("v_{2}^{a +}");
      hvncorr_m[i][j]->SetYTitle("v_{2}^{b -}");

      hangcorr[i][j] = new TH2D(Form("hangcorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorr[i][j]->SetDirectory(0);
      hangcorr[i][j]->Sumw2();
      hangcorr[i][j]->SetOption("colz");
      hangcorr[i][j]->SetXTitle("#Psi^{b -}");
      hangcorr[i][j]->SetYTitle("#Psi^{a +}");
      hangcorr_p[i][j] = new TH2D(Form("hangcorr_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorr_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorr_p[i][j]->SetDirectory(0);
      hangcorr_p[i][j]->Sumw2();
      hangcorr_p[i][j]->SetOption("colz");
      hangcorr_p[i][j]->SetXTitle("#Psi^{b -}");
      hangcorr_p[i][j]->SetYTitle("#Psi^{a +}");
      hangcorr_m[i][j] = new TH2D(Form("hangcorr_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorr_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorr_m[i][j]->SetDirectory(0);
      hangcorr_m[i][j]->Sumw2();
      hangcorr_m[i][j]->SetOption("colz");
      hangcorr_m[i][j]->SetXTitle("#Psi^{b -}");
      hangcorr_m[i][j]->SetYTitle("#Psi^{a +}");


      hangcorrHFpA[i][j] = new TH2D(Form("hangcorrHFpA_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpA_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorrHFpA[i][j]->SetDirectory(0);
      hangcorrHFpA[i][j]->Sumw2();
      hangcorrHFpA[i][j]->SetOption("colz");
      hangcorrHFpA[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpA[i][j]->SetYTitle("#Psi^{a +}");
      hangcorrHFpA_p[i][j] = new TH2D(Form("hangcorrHFpA_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpA_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorrHFpA_p[i][j]->SetDirectory(0);
      hangcorrHFpA_p[i][j]->Sumw2();
      hangcorrHFpA_p[i][j]->SetOption("colz");
      hangcorrHFpA_p[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpA_p[i][j]->SetYTitle("#Psi^{a +}");
      hangcorrHFpA_m[i][j] = new TH2D(Form("hangcorrHFpA_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpA_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorrHFpA_m[i][j]->SetDirectory(0);
      hangcorrHFpA_m[i][j]->Sumw2();
      hangcorrHFpA_m[i][j]->SetOption("colz");
      hangcorrHFpA_m[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpA_m[i][j]->SetYTitle("#Psi^{a +}");

      hangcorrHFpB[i][j] = new TH2D(Form("hangcorrHFpB_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpB_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorrHFpB[i][j]->SetDirectory(0);
      hangcorrHFpB[i][j]->Sumw2();
      hangcorrHFpB[i][j]->SetOption("colz");
      hangcorrHFpB[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpB[i][j]->SetYTitle("#Psi^{b -}");
      hangcorrHFpB_p[i][j] = new TH2D(Form("hangcorrHFpB_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpB_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorrHFpB_p[i][j]->SetDirectory(0);
      hangcorrHFpB_p[i][j]->Sumw2();
      hangcorrHFpB_p[i][j]->SetOption("colz");
      hangcorrHFpB_p[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpB_p[i][j]->SetYTitle("#Psi^{b -}");
      hangcorrHFpB_m[i][j] = new TH2D(Form("hangcorrHFpB_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpB_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorrHFpB_m[i][j]->SetDirectory(0);
      hangcorrHFpB_m[i][j]->Sumw2();
      hangcorrHFpB_m[i][j]->SetOption("colz");
      hangcorrHFpB_m[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpB_m[i][j]->SetYTitle("#Psi^{b -}");


      diffX[i][j] = new TH1D(Form("diffX_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffX_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),400,-1.0,1.0);
      diffX[i][j]->SetDirectory(0);
      diffX[i][j]->Sumw2();
      diffX[i][j]->SetOption("colz");
      diffX[i][j]->SetXTitle(Form("(v_{%d,x}^{obs,a}-v_{%d,x}^{obs,b})/2",Order,Order));
      diffX[i][j]->SetYTitle("Counts");
      diffX_p[i][j] = new TH1D(Form("diffX_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffX_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),400,-1.0,1.0);
      diffX_p[i][j]->SetDirectory(0);
      diffX_p[i][j]->Sumw2();
      diffX_p[i][j]->SetOption("colz");
      diffX_p[i][j]->SetXTitle(Form("(v_{%d,x}^{obs,a}-v_{%d,x}^{obs,b})/2",Order,Order));
      diffX_p[i][j]->SetYTitle("Counts");
      diffX_m[i][j] = new TH1D(Form("diffX_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffX_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),400,-1.0,1.0);
      diffX_m[i][j]->SetDirectory(0);
      diffX_m[i][j]->Sumw2();
      diffX_m[i][j]->SetOption("colz");
      diffX_m[i][j]->SetXTitle(Form("(v_{%d,x}^{obs,a}-v_{%d,x}^{obs,b})/2",Order,Order));
      diffX_m[i][j]->SetYTitle("Counts");

      diffY[i][j] = new TH1D(Form("diffY_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffY_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),400,-1.0,1.0);
      diffY[i][j]->SetDirectory(0);
      diffY[i][j]->Sumw2();
      diffY[i][j]->SetOption("colz");
      diffY[i][j]->SetXTitle(Form("(v_{%d,y}^{obs,a}-v_{%d,y}^{obs,b})/2",Order,Order));
      diffY[i][j]->SetYTitle("Counts");
      diffY_p[i][j] = new TH1D(Form("diffY_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffY_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),400,-1.0,1.0);
      diffY_p[i][j]->SetDirectory(0);
      diffY_p[i][j]->Sumw2();
      diffY_p[i][j]->SetOption("colz");
      diffY_p[i][j]->SetXTitle(Form("(v_{%d,y}^{obs,a}-v_{%d,y}^{obs,b})/2",Order,Order));
      diffY_p[i][j]->SetYTitle("Counts");
      diffY_m[i][j] = new TH1D(Form("diffY_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffY_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),400,-1.0,1.0);
      diffY_m[i][j]->SetDirectory(0);
      diffY_m[i][j]->Sumw2();
      diffY_m[i][j]->SetOption("colz");
      diffY_m[i][j]->SetXTitle(Form("(v_{%d,y}^{obs,a}-v_{%d,y}^{obs,b})/2",Order,Order));
      diffY_m[i][j]->SetYTitle("Counts");


      runqdifx[i][j] = new TH1D(Form("runqdifx_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdifx_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdifx[i][j]->SetDirectory(0);
      runqdifx[i][j]->Sumw2();
      runqdifx_p[i][j] = new TH1D(Form("runqdifx_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdifx_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdifx_p[i][j]->SetDirectory(0);
      runqdifx_p[i][j]->Sumw2();
      runqdifx_m[i][j] = new TH1D(Form("runqdifx_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdifx_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdifx_m[i][j]->SetDirectory(0);
      runqdifx_m[i][j]->Sumw2();


      runqdify[i][j] = new TH1D(Form("runqdify_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdify_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdify[i][j]->SetDirectory(0);
      runqdify[i][j]->Sumw2();
      runqdify_p[i][j] = new TH1D(Form("runqdify_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdify_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdify_p[i][j]->SetDirectory(0);
      runqdify_p[i][j]->Sumw2();
      runqdify_m[i][j] = new TH1D(Form("runqdify_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdify_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdify_m[i][j]->SetDirectory(0);
      runqdify_m[i][j]->Sumw2();

      runqdifcnt[i][j] = new TH1D(Form("runqdifcnt_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdifcnt_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdifcnt[i][j]->SetDirectory(0);
      runqdifcnt[i][j]->Sumw2();
      runqdifcnt_p[i][j] = new TH1D(Form("runqdifcnt_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdifcnt_p_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdifcnt_p[i][j]->SetDirectory(0);
      runqdifcnt_p[i][j]->Sumw2();
      runqdifcnt_m[i][j] = new TH1D(Form("runqdifcnt_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				Form("runqdifcnt_m_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nruns,runlist);
      runqdifcnt_m[i][j]->SetDirectory(0);
      runqdifcnt_m[i][j]->Sumw2();
      
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
	double xdiffoff_p = frame->GetXdiff_p(mapp[k][j],runno);
	double xdiffoff_m = frame->GetXdiff_m(mapp[k][j],runno);
	double ydiffoff = frame->GetYdiff(mapp[k][j],runno);
	double ydiffoff_p = frame->GetYdiff_p(mapp[k][j],runno);
	double ydiffoff_m = frame->GetYdiff_m(mapp[k][j],runno);

	double vnxp = frame->GetVnxEvt(mapp[k][j]);
	double vnxp_p = frame->GetVnxEvt_p(mapp[k][j]);
	double vnxp_m = frame->GetVnxEvt_m(mapp[k][j]);
	double vnxm = frame->GetVnxEvt(mapm[k][j]);
	double vnxm_p = frame->GetVnxEvt_p(mapm[k][j]);
	double vnxm_m = frame->GetVnxEvt_m(mapm[k][j]);
	double vnx = (vnxp+vnxm)/2.;
	double vnx_p = (vnxp_p+vnxm_p)/2.;
	double vnx_m = (vnxp_m+vnxm_m)/2.;

	double vnyp = frame->GetVnyEvt(mapp[k][j]);
	double vnyp_p = frame->GetVnyEvt_p(mapp[k][j]);
	double vnyp_m = frame->GetVnyEvt_m(mapp[k][j]);
	double vnym = frame->GetVnyEvt(mapm[k][j]);
	double vnym_p = frame->GetVnyEvt_p(mapm[k][j]);
	double vnym_m = frame->GetVnyEvt_m(mapm[k][j]);
	double vny = (vnyp+vnym)/2.;
	double vny_p = (vnyp_p+vnym_p)/2.;
	double vny_m = (vnyp_m+vnym_m)/2.;

	double vn = sqrt(pow(vnx,2)+pow(vny,2));
	double vnp = sqrt(pow(vnxp,2)+pow(vnyp,2));
	double vnm = sqrt(pow(vnxm,2)+pow(vnym,2));
	if(vnp>-2 && vnm>-2) {
	  //hvn[k][j]->Fill(vn);
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
	
	double vn_p = sqrt(pow(vnx_p,2)+pow(vny_p,2));
	double vnp_p = sqrt(pow(vnxp_p,2)+pow(vnyp_p,2));
	double vnm_p = sqrt(pow(vnxm_p,2)+pow(vnym_p,2));
	if(vnp_p>-2 && vnm_p>-2) {
	  //hvn_p[k][j]->Fill(vn_p);
	  hvncorr_p[k][j]->Fill(vnm_p,vnp_p);
	  hangcorr_p[k][j]->Fill(frame->GetAng_p(mapm[k][j]),frame->GetAng_p(mapp[k][j]));
	  hangcorrHFpA_p[k][j]->Fill(frame->GetAngHFp(mapp[k][j]),frame->GetAng_p(mapp[k][j]));
	  hangcorrHFpB_p[k][j]->Fill(frame->GetAngHFp(mapp[k][j]),frame->GetAng_p(mapm[k][j]));
	  double xdiff_p = (vnxp_p-vnxm_p)/2. - xdiffoff_p;
	  double ydiff_p = (vnyp_p-vnym_p)/2. - ydiffoff_p;
	  hDiff_p[k][j]->Fill(xdiff_p,ydiff_p); 
	  double vnobsx_p = 0;
	  double vnobsy_p = 0;
	  double vnobs_p = 0;
	  frame->Smear_p(mapp[k][j],vnx_p,vny_p,vnobsx_p,vnobsy_p,vnobs_p);
	  diffX_p[k][j]->Fill(xdiff_p,vnobs_p);
	  diffY_p[k][j]->Fill(ydiff_p,vnobs_p);
	  hResp_p[k][j]->Fill(vnobs_p,vn_p);
	  runqdifx_p[k][j]->Fill(runno,xdiff_p);
	  runqdify_p[k][j]->Fill(runno,ydiff_p);
	  runqdifcnt_p[k][j]->Fill(runno);
	}
	
	double vn_m = sqrt(pow(vnx_m,2)+pow(vny_m,2));
	double vnp_m = sqrt(pow(vnxp_m,2)+pow(vnyp_m,2));
	double vnm_m = sqrt(pow(vnxm_m,2)+pow(vnym_m,2));
	if(vnp_p>-2 && vnm_p>-2) {
	  //hvn_m[k][j]->Fill(vn_m);
	  hvncorr_m[k][j]->Fill(vnm_m,vnp_m);
	  hangcorr_m[k][j]->Fill(frame->GetAng_m(mapm[k][j]),frame->GetAng_m(mapp[k][j]));
	  hangcorrHFpA_m[k][j]->Fill(frame->GetAngHFp(mapp[k][j]),frame->GetAng_m(mapp[k][j]));
	  hangcorrHFpB_m[k][j]->Fill(frame->GetAngHFp(mapp[k][j]),frame->GetAng_m(mapm[k][j]));
	  double xdiff_m = (vnxp_m-vnxm_m)/2. - xdiffoff_m;
	  double ydiff_m = (vnyp_m-vnym_m)/2. - ydiffoff_m;
	  hDiff_m[k][j]->Fill(xdiff_m,ydiff_m); 
	  double vnobsx_m = 0;
	  double vnobsy_m = 0;
	  double vnobs_m = 0;
	  frame->Smear_m(mapp[k][j],vnx_m,vny_m,vnobsx_m,vnobsy_m,vnobs_m);
	  diffX_m[k][j]->Fill(xdiff_m,vnobs_m);
	  diffY_m[k][j]->Fill(ydiff_m,vnobs_m);
	  hResp_m[k][j]->Fill(vnobs_m,vn_m);
	  runqdifx_m[k][j]->Fill(runno,xdiff_m);
	  runqdify_m[k][j]->Fill(runno,ydiff_m);
	  runqdifcnt_m[k][j]->Fill(runno);
	}
      }
    }
    ++count;
    ++partcount;
    if(partcount == 10000) {
      cout<<count<<endl;
      partcount = 0;
    }
  }
  if(nfiles++<maxFiles && keepRunning == 1 && frame->AddFile()) goto rep;
  //////////////
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
    ((TH2D *)frame->Get2d(i))->Write(Form("vn2d%s",frame->GetSide(i).data()));
    ((TH1D *)frame->Get1d(i))->Write(Form("vn1d%s",frame->GetSide(i).data()));
    TH1D * tmp = (TH1D *) frame->Get1dMult(i);
    tmp->Divide((TH1D *)frame->Get1d(i));
    tmp->Write(Form("vn1dMult%s",frame->GetSide(i).data()));
    ((TH1D *)frame->GetMult(i))->Write(Form("mult%s",frame->GetSide(i).data()));
    int cb = centloc[i];
    int pb = ptloc[i];
    if(hvn[cb][pb]==0) {
      hvn[cb][pb] = frame->Get1d(i);
    } else {
      hvn[cb][pb]->Add(frame->Get1d(i));
    }
    if(pb < 12) {
      double ptav = frame->GetAvPt(i);
      double vnv = frame->GetVn(i);
      double vnve = frame->GetVnErr(i);
      int bin = hvnpt_pos_eta[0]->FindBin(ptav);
      if(frame->GetSide(i).find("pos_eta")!=std::string::npos) {
	hvnpt_pos_eta[cb]->SetBinContent(bin,vnv);
	hvnpt_pos_eta[cb]->SetBinError(bin,vnve);
      } else if (frame->GetSide(i).find("neg_eta")!=std::string::npos){
	hvnpt_neg_eta[cb]->SetBinContent(bin,vnv);
	hvnpt_neg_eta[cb]->SetBinError(bin,vnve);
      }
    }
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
    subdirs[i]->cd();
    hvnpt_pos_eta[i]->Write("vnpt_pos_eta");
    hvnpt_neg_eta[i]->Write("vnpt_neg_eta");
  }
  fout->cd();
  //////////////
  TFile * fout_p = new TFile("checkXY_p.root","recreate");
  fout_p->cd();
  TDirectory * subdirs_p[ncentbins];
  TDirectory * subsubdirs_p[ncentbins][nptbins];
  for(int i = 0; i<ncentbins; i++) {
    subdirs_p[i] = fout_p->mkdir(Form("%d_%d",centbins[i],centbins[i+1]));
    subdirs_p[i]->cd();
    TH1D * sp_p = frame->GetSpectra(ROIp[0][i]);
    sp_p->Write();
    for(int j = 0; j<nptbins; j++) {
      subsubdirs_p[i][j] = subdirs_p[i]->mkdir(Form("%04.2f_%04.2f",ptbinsMin[j],ptbinsMax[j]));
    }
  }
  for(int i = 0; i<frame->GetNrange(); i++) {
    subsubdirs_p[centloc[i]][ptloc[i]]->cd();
    ((TH2D *)frame->Get2d_p(i))->Write(Form("vn2d%s",frame->GetSide(i).data()));
    ((TH1D *)frame->Get1d_p(i))->Write(Form("vn1d%s",frame->GetSide(i).data()));
    TH1D * tmp_p = (TH1D *) frame->Get1dMult_p(i);
    tmp_p->Divide((TH1D *)frame->Get1d_p(i));
    tmp_p->Write(Form("vn1dMult%s",frame->GetSide(i).data()));
    ((TH1D *)frame->GetMult_p(i))->Write(Form("mult%s",frame->GetSide(i).data()));
    int cb = centloc[i];
    int pb = ptloc[i];
     if(hvn_p[cb][pb]==0) {
      hvn_p[cb][pb] = frame->Get1d_p(i);
    } else {
      hvn_p[cb][pb]->Add(frame->Get1d_p(i));
    }
    if(pb < 12) {
      double ptav_p = frame->GetAvPt_p(i);
      double vnv_p = frame->GetVn_p(i);
      double vnve_p = frame->GetVnErr_p(i);
      int bin = hvnpt_pos_eta_p[0]->FindBin(ptav_p);

      if(frame->GetSide(i).find("pos_eta")!=std::string::npos) {
	hvnpt_pos_eta_p[cb]->SetBinContent(bin,vnv_p);
	hvnpt_pos_eta_p[cb]->SetBinError(bin,vnve_p);
      } else if (frame->GetSide(i).find("neg_eta")!=std::string::npos){
	hvnpt_neg_eta_p[cb]->SetBinContent(bin,vnv_p);
	hvnpt_neg_eta_p[cb]->SetBinError(bin,vnve_p);
      }
    }
  }
    
  for(int i = 0; i<ncentbins; i++) {
    for(int j = 0; j<nptbins; j++) {
      subsubdirs_p[i][j]->cd();
      hvn_p[i][j]->Write("vn");
      hvncorr_p[i][j]->Write("vncorr");
      hangcorr_p[i][j]->Write("angcorr");
      hangcorrHFpA_p[i][j]->Write("angcorr_HFpA");
      hangcorrHFpB_p[i][j]->Write("angcorr_HFpB");
      hDiff_p[i][j]->Write("diff2d");
      hResp_p[i][j]->Write("resp2d");
      diffX_p[i][j]->Write("diffX");
      diffY_p[i][j]->Write("diffY");
      runqdifx_p[i][j]->Write("qdifx");
      runqdify_p[i][j]->Write("qdify");
      runqdifcnt_p[i][j]->Write("qdifcnt");
      CreatePlots(i,j,"_p");
    }
    subdirs_p[i]->cd();
    hvnpt_pos_eta_p[i]->Write("vnpt_pos_eta_p");
    hvnpt_neg_eta_p[i]->Write("vnpt_neg_eta_p");
  }
  fout_p->cd();
  //////////////
  TFile * fout_m = new TFile("checkXY_m.root","recreate");
  fout_m->cd();
  TDirectory * subdirs_m[ncentbins];
  TDirectory * subsubdirs_m[ncentbins][nptbins];
  for(int i = 0; i<ncentbins; i++) {
    subdirs_m[i] = fout_m->mkdir(Form("%d_%d",centbins[i],centbins[i+1]));
    subdirs_m[i]->cd();
    TH1D * sp_m = frame->GetSpectra(ROIp[0][i]);
    sp_m->Write();
    for(int j = 0; j<nptbins; j++) {
      subsubdirs_m[i][j] = subdirs_m[i]->mkdir(Form("%04.2f_%04.2f",ptbinsMin[j],ptbinsMax[j]));
    }
  }
  for(int i = 0; i<frame->GetNrange(); i++) {
    subsubdirs_m[centloc[i]][ptloc[i]]->cd();
    ((TH2D *)frame->Get2d_m(i))->Write(Form("vn2d%s",frame->GetSide(i).data()));
    ((TH1D *)frame->Get1d_m(i))->Write(Form("vn1d%s",frame->GetSide(i).data()));
    TH1D * tmp_m = (TH1D *) frame->Get1dMult_m(i);
    tmp_m->Divide((TH1D *)frame->Get1d_m(i));
    tmp_m->Write(Form("vn1dMult%s",frame->GetSide(i).data()));
    ((TH1D *)frame->GetMult_m(i))->Write(Form("mult%s",frame->GetSide(i).data()));
    int cb = centloc[i];
    int pb = ptloc[i];
     if(hvn_m[cb][pb]==0) {
      hvn_m[cb][pb] = frame->Get1d_m(i);
    } else {
      hvn_m[cb][pb]->Add(frame->Get1d_m(i));
    }
    if(pb < 12) {
      double ptav_m = frame->GetAvPt_m(i);
      double vnv_m = frame->GetVn_m(i);
      double vnve_m = frame->GetVnErr_m(i);
      int bin = hvnpt_pos_eta_m[0]->FindBin(ptav_m);

      if(frame->GetSide(i).find("pos_eta")!=std::string::npos) {
	hvnpt_pos_eta_m[cb]->SetBinContent(bin,vnv_m);
	hvnpt_pos_eta_m[cb]->SetBinError(bin,vnve_m);
      } else if (frame->GetSide(i).find("neg_eta")!=std::string::npos){
	hvnpt_neg_eta_m[cb]->SetBinContent(bin,vnv_m);
	hvnpt_neg_eta_m[cb]->SetBinError(bin,vnve_m);
      }
    }
  }
    
  for(int i = 0; i<ncentbins; i++) {
    for(int j = 0; j<nptbins; j++) {
      subsubdirs_m[i][j]->cd();
      hvn_m[i][j]->Write("vn");
      hvncorr_m[i][j]->Write("vncorr");
      hangcorr_m[i][j]->Write("angcorr");
      hangcorrHFpA_m[i][j]->Write("angcorr_HFpA");
      hangcorrHFpB_m[i][j]->Write("angcorr_HFpB");
      hDiff_m[i][j]->Write("diff2d");
      hResp_m[i][j]->Write("resp2d");
      diffX_m[i][j]->Write("diffX");
      diffY_m[i][j]->Write("diffY");
      runqdifx_m[i][j]->Write("qdifx");
      runqdify_m[i][j]->Write("qdify");
      runqdifcnt_m[i][j]->Write("qdifcnt");
      CreatePlots(i,j,"_m");
    }
    subdirs_m[i]->cd();
    hvnpt_pos_eta_m[i]->Write("vnpt_pos_eta_m");
    hvnpt_neg_eta_m[i]->Write("vnpt_neg_eta_m");
  }
  fout_m->cd();

  for(int i = 0; i<frame->GetNrange(); i++) {
    frame->SaveFrameworkRuns(i,subsubdirs[centloc[i]][ptloc[i]]);
    frame->SaveFrameworkRuns_p(i,subsubdirs_p[centloc[i]][ptloc[i]]);
    frame->SaveFrameworkRuns_m(i,subsubdirs_m[centloc[i]][ptloc[i]]);
  }

}

void   CreatePlots(int icent, int ipt, string suffix){

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
  
  clow->Print(Form("figures/checkXY/%s%s.pdf",clow->GetName(),suffix.data()),"pdf");
  clow->Close();
  gSystem->ProcessEvents();
  delete clow;
  clow = 0;
  return;
}
