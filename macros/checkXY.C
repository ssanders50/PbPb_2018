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
#include <ctime>
#include "../CMSSW_10_3_3_patch1/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;
using namespace std;
static const int maxFiles = 1000;
static const int maxEvents = 200000000;
static volatile int keepRunning = 1;
void intHandler(int dummy){
  keepRunning = 0;
}
int CS = 0;
TFile * tf=NULL;
static const int nptbins = 15;
static const int nhistbins = 100;
static const double maxvn = 0.6;
static const float ptbinsMin[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0, 2.5,  3.0,  3.5, 0.3, 0.6, 1.0};
static const float ptbinsMax[]={0.4, 0.5, 0.6,  0.8,  1.0, 1.25,  1.50,   2.0,  2.5, 3.0,  3.5,  4.0, 3.0, 3.0, 3.0};
static const int ncentbins = 13;
static const int centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,100};
float sub1[2] = {-1.2,0};
float sub2[2] = {0,1.2};
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
TH1D * hvn[ncentbins][nptbins]={0};
TH1D * diffX[ncentbins][nptbins];
TH1D * diffY[ncentbins][nptbins];
TH1D * runqdifx[ncentbins][nptbins];
TH1D * runqdify[ncentbins][nptbins];
TH1D * runqdifcnt[ncentbins][nptbins];
TH1D * hvnpt_sub1[ncentbins];
TH1D * hvnpt_sub2[ncentbins];
TH1D * hvnpt_full[ncentbins];
TH1D * hCent[ncentbins][nptbins];
TH1D * hVtx[ncentbins][nptbins];
TH1I * hruns;
double runlist[500];
int nruns = 0;
int ROI[ncentbins][nptbins];
int maproi[ncentbins][nptbins];
int centloc[500];
int ptloc[500];
Framework * frame;

void CreatePlots(int icent, int ipt, string suffix = "");

void checkXY(string inlist="filelist.dat",int cs = 0){
  cout<<"==========================================="<<endl;
  cout<<"Kill detached process with: kill -2 [id]"<<endl;
  cout<<"==========================================="<<endl;
  CS = cs;
  std::clock_t start = 0;
  std::clock_t last = 0;
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
  frame = new Framework(inlist,true, nhistbins, maxvn, CS);
  frame->SetMinMult(2);
  frame->SetRuns(nruns,runlist);
  if(CS==0) {
    if(!frame->LoadOffsets("checkXY_offsets.root")) {
      cout<<"Failed to open offsets file"<<endl; }
  } else if (CS==1) {
    if(!frame->LoadOffsets("checkXY_p_offsets.root")) {
      cout<<"Failed to open offsets_p file"<<endl; }
  } else if (CS==-1) {
    if(!frame->LoadOffsets("checkXY_m_offsets.root")) {
      cout<<"Failed to open offsets_m file"<<endl; }
  }
  int iloc = 0;
  for(int i = 0; i<ncentbins; i++) {
    for(int j = 0; j<nptbins; j++) {
      ROI[i][j]=frame->SetROIRange(Order,centbins[i],centbins[i+1],sub1,sub2, ptbinsMin[j], ptbinsMax[j]);
      frame->SetROIEP(ROI[i][j],HFm2,HFp2,trackmid2);
      maproi[i][j]=iloc;
      ptloc[iloc] = j;
      centloc[iloc++]=i;

      hCent[i][j] = new TH1D(Form("cent_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("cent_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),200,0,100);
      hCent[i][j]->SetDirectory(0);
      hCent[i][j]->Sumw2();
      hCent[i][j]->SetOption("colz");
      hCent[i][j]->SetXTitle("Centrality");
      hCent[i][j]->SetYTitle("Count");

      hVtx[i][j] = new TH1D(Form("vtx_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("vtx_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),200,-25,25);
      hVtx[i][j]->SetDirectory(0);
      hVtx[i][j]->Sumw2();
      hVtx[i][j]->SetOption("colz");
      hVtx[i][j]->SetXTitle("Vertex");
      hVtx[i][j]->SetYTitle("Count");


      hDiff[i][j] = new TH2D(Form("diff_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("diff_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,-1,1,nhistbins,-1,1);
      hDiff[i][j]->SetDirectory(0);
      hDiff[i][j]->Sumw2();
      hDiff[i][j]->SetOption("colz");
      hDiff[i][j]->SetXTitle(Form("(v_{%d,x}^{obs,a}-v_{%d,x}^{obs,b})/2",Order,Order));
      hDiff[i][j]->SetYTitle(Form("(v_{%d,y}^{obs,a}-v_{%d,y}^{obs,b})/2",Order,Order));
   

      hResp[i][j] = new TH2D(Form("resp_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			     Form("resp_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn,nhistbins,0,maxvn);
      hResp[i][j]->SetDirectory(0);
      hResp[i][j]->Sumw2();
      hResp[i][j]->SetOption("colz");
      hResp[i][j]->SetXTitle(Form("v_{%d}^{obs}",Order));
      hResp[i][j]->SetYTitle(Form("v_{%d}^{obs}",Order));

      hvnpt_sub1[i] = new TH1D(Form("hvnpt_sub1_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_sub1_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_sub1[i]->SetDirectory(0);
      hvnpt_sub1[i]->Sumw2();
      hvnpt_sub1[i]->SetOption("colz");
      hvnpt_sub1[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_sub1[i]->SetYTitle(Form("v_{%d} (%03.1f#le #eta < %03.1f)",Order,sub1[0],sub1[1]));


      hvnpt_sub2[i] = new TH1D(Form("hvnpt_sub2_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_sub2_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_sub2[i]->SetDirectory(0);
      hvnpt_sub2[i]->Sumw2();
      hvnpt_sub2[i]->SetOption("colz");
      hvnpt_sub2[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_sub2[i]->SetYTitle(Form("v_{%d} (%03.1f#le #eta < %03.1f)",Order,sub2[0],sub2[1]));

      hvnpt_full[i] = new TH1D(Form("hvnpt_full_%d_%d",centbins[i],centbins[i+1]),
			   Form("hvnpt_full_%d_%d",centbins[i],centbins[i+1]),500,0,5);
      hvnpt_full[i]->SetDirectory(0);
      hvnpt_full[i]->Sumw2();
      hvnpt_full[i]->SetOption("colz");
      hvnpt_full[i]->SetXTitle("p_{T} (GeV/c)");
      hvnpt_full[i]->SetYTitle(Form("v_{%d}",Order));


      hvncorr[i][j] = new TH2D(Form("hvncorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hvncorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),nhistbins,0,maxvn,nhistbins,0,maxvn);
      hvncorr[i][j]->SetDirectory(0);
      hvncorr[i][j]->Sumw2();
      hvncorr[i][j]->SetOption("colz");
      hvncorr[i][j]->SetXTitle("v_{2}^{a +}");
      hvncorr[i][j]->SetYTitle("v_{2}^{b -}");

      hangcorr[i][j] = new TH2D(Form("hangcorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorr_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorr[i][j]->SetDirectory(0);
      hangcorr[i][j]->Sumw2();
      hangcorr[i][j]->SetOption("colz");
      hangcorr[i][j]->SetXTitle("#Psi^{b -}");
      hangcorr[i][j]->SetYTitle("#Psi^{a +}");


      hangcorrHFpA[i][j] = new TH2D(Form("hangcorrHFpA_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpA_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorrHFpA[i][j]->SetDirectory(0);
      hangcorrHFpA[i][j]->Sumw2();
      hangcorrHFpA[i][j]->SetOption("colz");
      hangcorrHFpA[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpA[i][j]->SetYTitle("#Psi^{a +}");

      hangcorrHFpB[i][j] = new TH2D(Form("hangcorrHFpB_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			       Form("hangcorrHFpB_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
				nhistbins,-2,2,nhistbins,-2,2);
      hangcorrHFpB[i][j]->SetDirectory(0);
      hangcorrHFpB[i][j]->Sumw2();
      hangcorrHFpB[i][j]->SetOption("colz");
      hangcorrHFpB[i][j]->SetXTitle("#Psi^{HF+}");
      hangcorrHFpB[i][j]->SetYTitle("#Psi^{b -}");

      diffX[i][j] = new TH1D(Form("diffX_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffX_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),800,-1.0,1.0);
      diffX[i][j]->SetDirectory(0);
      diffX[i][j]->Sumw2();
      diffX[i][j]->SetOption("colz");
      diffX[i][j]->SetXTitle(Form("(v_{%d,x}^{obs,a}-v_{%d,x}^{obs,b})/2",Order,Order));
      diffX[i][j]->SetYTitle("Counts");

      diffY[i][j] = new TH1D(Form("diffY_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),
			      Form("diffY_%d_%d_%04.2f_%04.2f",centbins[i],centbins[i+1],ptbinsMin[j],ptbinsMax[j]),800,-1.0,1.0);
      diffY[i][j]->SetDirectory(0);
      diffY[i][j]->Sumw2();
      diffY[i][j]->SetOption("colz");
      diffY[i][j]->SetXTitle(Form("(v_{%d,y}^{obs,a}-v_{%d,y}^{obs,b})/2",Order,Order));
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
  start = std::clock();
 rep: for(int i = 0; i<frame->GetN(); i++) {
    signal(SIGINT,intHandler);
    if(keepRunning == 0) break;
    int runno = frame->GetRunno(i);
    int irun = hruns->FindBin(runno+0.1);
    frame->AddEvent(i);
    for(int k = 0; k<ncentbins; k++) {
      for(int j = 0; j<nptbins; j++) {
	if(frame->GetVnxEvt_sub1(maproi[k][j])<-1.) continue;
	if(frame->GetVnxEvt_sub2(maproi[k][j])<-1.) continue;
	if(frame->GetVnyEvt_sub1(maproi[k][j])<-1.) continue;
	if(frame->GetVnyEvt_sub2(maproi[k][j])<-1.) continue;
	double xdiffoff_sub1 = frame->GetXdiff_sub1(maproi[k][j],runno);
	double ydiffoff_sub1 = frame->GetYdiff_sub2(maproi[k][j],runno);
	
	double vnx_sub1 = frame->GetVnxEvt_sub1(maproi[k][j]);
	double vnx_sub2 = frame->GetVnxEvt_sub2(maproi[k][j]);
	double vnx_full =  frame->GetVnxEvt_full(maproi[k][j]);
	
	double vny_sub1 = frame->GetVnyEvt_sub1(maproi[k][j]);
	double vny_sub2 = frame->GetVnyEvt_sub2(maproi[k][j]);
	double vny_full =  frame->GetVnyEvt_full(maproi[k][j]);
	
	double vn_sub1 = sqrt(pow(vnx_sub1,2)+pow(vny_sub1,2));
	double vn_sub2 = sqrt(pow(vnx_sub2,2)+pow(vny_sub2,2));
	double vn_full = sqrt(pow(vnx_full,2)+pow(vny_full,2));
	if(vn_sub1>-2 && vn_sub2>-2) {
	  hvncorr[k][j]->Fill(vn_sub1,vn_sub2);
	  hangcorr[k][j]->Fill(frame->GetAng_sub1(maproi[k][j]),frame->GetAng_sub2(maproi[k][j]));
	  // 	  hangcorrHFpA[k][j]->Fill(frame->GetAngHFp(mapp[k][j]),frame->GetAng(mapp[k][j]));
	  // 	  hangcorrHFpB[k][j]->Fill(frame->GetAngHFp(mapp[k][j]),frame->GetAng(mapm[k][j]));
	// 	  double xdiff = (vnxp-vnxm)/2. - xdiffoff;
	// 	  double ydiff = (vnyp-vnym)/2. - ydiffoff;
	// 	  hDiff[k][j]->Fill(xdiff,ydiff); 
	// 	  double vnobsx = 0;
	// 	  double vnobsy = 0;
	// 	  double vnobs = 0;
	// 	  frame->Smear(mapp[k][j],vnx,vny,vnobsx,vnobsy,vnobs);
	// 	  diffX[k][j]->Fill(xdiff,vnobs);
	// 	  diffY[k][j]->Fill(ydiff,vnobs);
	// 	  hResp[k][j]->Fill(vnobs,vn);
	// 	  runqdifx[k][j]->Fill(runno,xdiff);
	// 	  runqdify[k][j]->Fill(runno,ydiff);
	// 	  runqdifcnt[k][j]->Fill(runno);
	  hCent[k][j]->Fill(frame->GetCent(maproi[k][j]));
	  hVtx[k][j]->Fill(frame->GetVtx(maproi[k][j]));
	}
	
      }
    }
    ++count;
    ++partcount;
    if(partcount == 100000) {
      std::clock_t now = std::clock();
      double tpe = (now - last) / (double)CLOCKS_PER_SEC / 50.;
      last = now;
      double duration = (now - start) / (double)CLOCKS_PER_SEC;
      duration/=60.;
      double hr = duration/60.;
      cout<<setw(10)<<count<<"   time per event (ms): "<<setw(15)<<setprecision(5)<<fixed<<tpe<<"   elapsed time (min/h): "<<setw(15)<<setprecision(5)<<fixed<<duration<<setw(15)<<setprecision(5)<<fixed<<hr<<endl;
      partcount = 0;
    }
    if(count>maxEvents) break;
  }
  if(nfiles++<maxFiles && keepRunning == 1 && frame->AddFile() && count<maxEvents) goto rep;
 //  //////////////
  TFile * fout = NULL;
  if(CS == 0) {
    fout = new TFile("checkXY.root","recreate");
  } else if (CS==1) {
    fout = new TFile("checkXY_p.root","recreate");
  } else if (CS==-1) {
    fout = new TFile("checkXY_m.root","recreate");
  }
  fout->cd();
  TDirectory * subdirs[ncentbins];
  TDirectory * subsubdirs[ncentbins][nptbins];
  for(int i = 0; i<ncentbins; i++) {
    subdirs[i] = fout->mkdir(Form("%d_%d",centbins[i],centbins[i+1]));
    subdirs[i]->cd();
    TH1D * sp = frame->GetSpectra(ROI[0][i]);
    sp->Write();
    for(int j = 0; j<nptbins; j++) {
      subsubdirs[i][j] = subdirs[i]->mkdir(Form("%04.2f_%04.2f",ptbinsMin[j],ptbinsMax[j]));
    }
  }
  for(int i = 0; i<frame->GetNrange(); i++) {
    subsubdirs[centloc[i]][ptloc[i]]->cd();
    ((TH2D *)frame->Get2d_sub1(i))->Write("vn2d_sub1");
    ((TH1D *)frame->Get1d_sub1(i))->Write("vn1d_sub1");
    ((TH2D *)frame->Get2d_sub2(i))->Write("vn2d_sub2");
    ((TH1D *)frame->Get1d_sub2(i))->Write("vn1d_sub2");
    ((TH2D *)frame->Get2d_full(i))->Write("vn2d_full");
    ((TH1D *)frame->Get1d_full(i))->Write("vn1d_full");
    TH1D * tmp = (TH1D *) frame->Get1dMult_sub1(i);
    tmp->Divide((TH1D *)frame->Get1d_sub1(i));
    tmp->Write("vn1dMult_sub1");
    ((TH1D *)frame->GetMult_sub1(i))->Write("mult_sub1");

    tmp = (TH1D *) frame->Get1dMult_sub2(i);
    tmp->Divide((TH1D *)frame->Get1d_sub2(i));
    tmp->Write("vn1dMult_sub2");
    ((TH1D *)frame->GetMult_sub2(i))->Write("mult_sub2");

    tmp = (TH1D *) frame->Get1dMult_full(i);
    tmp->Divide((TH1D *)frame->Get1d_full(i));
    tmp->Write("vn1dMult_full");
    ((TH1D *)frame->GetMult_full(i))->Write("mult_full");

 //    int cb = centloc[i];
 //    int pb = ptloc[i];
 //    if(hvn[cb][pb]==0) {
 //      hvn[cb][pb] = frame->Get1d(i);
 //    } else {
 //      hvn[cb][pb]->Add(frame->Get1d(i));
 //    }
 //    if(pb < 12) {
 //      double ptav = frame->GetAvPt(i);
 //      double vnv = frame->GetVn(i);
 //      double vnve = frame->GetVnErr(i);
 //      int bin = hvnpt_pos_eta[0]->FindBin(ptav);
 //      if(frame->GetSide(i).find("pos_eta")!=std::string::npos) {
 // 	hvnpt_pos_eta[cb]->SetBinContent(bin,vnv);
 // 	hvnpt_pos_eta[cb]->SetBinError(bin,vnve);
 //      } else if (frame->GetSide(i).find("neg_eta")!=std::string::npos){
 // 	hvnpt_neg_eta[cb]->SetBinContent(bin,vnv);
 // 	hvnpt_neg_eta[cb]->SetBinError(bin,vnve);
 //      }
 //    }
 //  }
 //  for(int i = 0; i<ncentbins; i++) {
 //    for(int j = 0; j<nptbins; j++) {
 //      subsubdirs[i][j]->cd();
 //      hvn[i][j]->Write("vn");
 //      hvncorr[i][j]->Write("vncorr");
 //      hangcorr[i][j]->Write("angcorr");
 //      hangcorrHFpA[i][j]->Write("angcorr_HFpA");
 //      hangcorrHFpB[i][j]->Write("angcorr_HFpB");
 //      hDiff[i][j]->Write("diff2d");
 //      hResp[i][j]->Write("resp2d");
 //      diffX[i][j]->Write("diffX");
 //      diffY[i][j]->Write("diffY");
 //      runqdifx[i][j]->Write("qdifx");
 //      runqdify[i][j]->Write("qdify");
 //      runqdifcnt[i][j]->Write("qdifcnt");
 //      hCent[i][j]->Write("cent");
 //      hVtx[i][j]->Write("vtx");
 //      CreatePlots(i,j);
 //    }
 //    subdirs[i]->cd();
 //    hvnpt_pos_eta[i]->Write("vnpt_pos_eta");
 //    hvnpt_neg_eta[i]->Write("vnpt_neg_eta");
  }
 //  fout->cd();

 //  for(int i = 0; i<frame->GetNrange(); i++) {
 //    frame->SaveFrameworkRuns(i,subsubdirs[centloc[i]][ptloc[i]]);
 //  }

}

// void   CreatePlots(int icent, int ipt, string suffix){

//   string cname = Form("xy_%d_%d_%3.1f_%3.1f",centbins[icent],centbins[icent+1],ptbinsMin[ipt],ptbinsMax[ipt]);
//   TCanvas * clow = new TCanvas(cname.data(),cname.data(),1300,1300);
//   clow->Divide(2,2);
//   clow->cd(1);

//   frame->Get2d(mapm[icent][ipt])->Draw();
//   gPad->SetGrid(1,1);
//   TLine * l1a = new TLine(-1.4,0,1.4,0);
//   l1a->SetLineWidth(2);
//   l1a->Draw();
//   TLine * l2a = new TLine(0,-1.4,0,1.4);
//   l2a->SetLineWidth(2);
//   l2a->Draw();
//   TEllipse * el1a = new TEllipse(0,0,1,1);
//   el1a->SetFillStyle(0);
//   el1a->SetLineWidth(3);
//   el1a->SetLineStyle(2);
//   el1a->SetLineColor(kRed);
//   el1a->Draw();
//   clow->cd(2);
//   frame->Get2d(mapp[icent][ipt])->Draw();
//   gPad->SetGrid(1,1);
//   TLine * l1b = new TLine(-1.4,0,1.4,0);
//   l1b->SetLineWidth(2);
//   l1b->Draw();
//   TLine * l2b = new TLine(0,-1.4,0,1.4);
//   l2b->SetLineWidth(2);
//   l2b->Draw();
//   TEllipse * el1b = new TEllipse(0,0,1,1);
//   el1b->SetFillStyle(0);
//   el1b->SetLineWidth(3);
//   el1b->SetLineStyle(2);
//   el1b->SetLineColor(kRed);
//   el1b->Draw();
//   clow->cd(3);
//   hDiff[icent][ipt]->Draw();
//   gPad->SetGrid(1,1);
//   TLine * l1c = new TLine(-1.4,0,1.4,0);
//   l1c->SetLineWidth(2);
//   l1c->Draw();
//   TLine * l2c = new TLine(0,-1.4,0,1.4);
//   l2c->SetLineWidth(2);
//   l2c->Draw();
//   clow->cd(4);
//   TPaveText * tp = new TPaveText(0.2,0.2,0.8,0.8);
//   tp->SetBorderSize(0);
//   tp->SetFillColor(kWhite);
//   tp->AddText("2018 PbPb");
//   tp->AddText(Form("%d - %d%c",centbins[icent],centbins[icent+1],'%'));
//   tp->AddText(Form("%3.1f < p_{T} <  %3.1f GeV/c",ptbinsMin[ipt],ptbinsMax[ipt]));
//   tp->Draw();
  
//   clow->Print(Form("figures/checkXY/%s%s.pdf",clow->GetName(),suffix.data()),"pdf");
//   clow->Close();
//   gSystem->ProcessEvents();
//   delete clow;
//   clow = 0;
//   return;
// }
