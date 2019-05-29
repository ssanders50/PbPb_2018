#ifndef FRAMEWORK
#define FRAMEWORK
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
static const uint minRunRange = 326381;
static const uint maxRunRange = 327565;
static const int MAXROI = 500;
static const int MAXRUNS = 500; 
class Framework{
public:
  Framework(string filelist="filelist.dat");
  bool AddFile();
  int GetN(){return maxevents;}
  void GetEntry(int i){tr->GetEntry(i);}
  int GetNrange(){return nrange;}
  int SetROIRange(int order, int minCent, int maxCent, double minEta, double maxEta, double minPt, double maxPt);
  void ShowAllROIRanges();
  void SetROIEP(int roi, int EPA, int EPB, int EPC, int EPA2=-1, int EPB2=-1, int EPC2 = -1);
  void AddEvent(int evt);
  double GetAvPt(int roi){return r[roi].pt/r[roi].wn;}
  double GetqnA(int roi){return r[roi].qn/r[roi].wnA;}
  double GetVn(int roi) {return GetqnA(roi)/GetqABC(roi);}
  double GetVnErrSubEvt(int roi);
  double GetVnErr(int roi) {return GetqnAError(roi)/GetqABC(roi);}
  double GetqABC(int roi);
  TH1D * GetSpectra(int roi, bool effCorrect = true);
  TH2D * Get2d(int roi){return r[roi].vn2d;}
  TH1D * Get1d(int roi){return r[roi].vn1d;}
  TH1D * Get1dMult(int roi){return r[roi].vn1dMult;}
  int GetMinCent(int roi){return r[roi].minCent;}
  int GetMaxCent(int roi){return r[roi].maxCent;}
  int GetMinMult(){return minMult;}
  int GetMaxMult(){return maxMult;}
  int GetMinRun(){return minRun;}
  int GetMaxRun(){return maxRun;}

  TH1D * GetRuns(){return runs;}
  TH1D * GetMult(int roi){return r[roi].mult;}
  double GetVnxEvt(int roi){return vnxEvt[roi];}
  double GetVnyEvt(int roi){return vnyEvt[roi];}
  uint GetRunno(int evt){GetEntry(evt); return runno_;}
  void SetMinMult(int val){minMult = val;}
  void SetMaxMult(int val){maxMult = val;}
  void SetMinRun(int val){minRun = val;}
  void SetMaxRun(int val){maxRun = val;}
  void SetRuns(int nruns, double * runs); 
  void SaveFrameworkRuns(int roi, TDirectory * dir);
  bool LoadOffsets(string offname);
  double GetXdiff(int roi, int runno);
  double GetYdiff(int roi, int runno);
  FILE * flist;
private:
  TH1D * runbins=NULL;
  TRandom * ran;
  double vnxEvt[MAXROI];
  double vnyEvt[MAXROI];
  double GetVnSub(int roi,int i) {return r[roi].qnSub[i]/r[roi].wnASub[i]/GetqABC(roi);}
  double GetqnAError(int roi){return sqrt(r[roi].qne)/r[roi].wnA;}
  int maxevents;
  TFile * tf;
  TFile *foff;
  TTree * tr;
  double centval;
  int ntrkval;
  int Noff;
  int minMult = -1;
  int maxMult = -1;
  uint minRun = 0;
  uint maxRun = 0;
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
  TH1D * runs=NULL;
  TH1D * runcnt;
  TH1D * runqx;
  TH1D * runqy;
  TH1D * qdifx;
  TH1D * qdify;
  int nruns;
  double runlist[MAXRUNS];

  struct range {
    int order;
    int orderIndx;
    int minCent;
    int maxCent;
    double minEta;
    double maxEta;
    double minPt;
    double maxPt;
    int minEtaBin;
    int maxEtaBin;
    int minPtBin;
    int maxPtBin;
    int A=-1;
    int B=-1;
    int C=-1;
    int A2=-1;
    int B2=-1;
    int C2=-1;
    double qn=0;
    double qne=0;
    double pt=0;
    double wn=0;
    double wne=0;
    double wnA = 0;
    double qAB = 0;
    double wAB = 0;
    double qAC = 0;
    double wAC = 0;
    double qBC = 0;
    double wBC = 0;
    double qnSub[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSube[10]={0,0,0,0,0,0,0,0,0,0};
    double wnASub[10]={0,0,0,0,0,0,0,0,0,0};
    TH2D * vn2d;
    TH1D * vn1d;
    TH1D * vn1dMult;
    TH1D * mult;
    TH1D * runcnt;
    TH1D * runqx;
    TH1D * runqy;
    TH1D * qoffx;
    TH1D * qoffy;
    TH1D * qdifx;
    TH1D * qdify;
  } r[500];
  int nrange = 0;
};
bool Framework::LoadOffsets(string offname){
  foff = new TFile(offname.data(),"read"); 
  if(foff->IsZombie()) return false;

  return true;
}
void Framework::SetRuns(int nr, double * runl){
  runs = new TH1D("FrameworkRuns","FrameworkRuns",nr,runl);
  runs->SetDirectory(0);
  nruns = nr;
  for(int i = 0; i<=nr; i++) runlist[i]=runl[i];
}

  double Framework::GetXdiff(int roi, int runno){
    if(foff->IsZombie()) return 0;
    int runbin = runs->GetBin(runno);
    return r[roi].qdifx->GetBinContent(runbin);
  }
  double Framework::GetYdiff(int roi, int runno){
    if(foff->IsZombie()) return 0;
    int runbin = runs->GetBin(runno);
    return r[roi].qdify->GetBinContent(runbin);
  }
void Framework::SaveFrameworkRuns(int roi, TDirectory * dir){
  TDirectory * savedir = gDirectory;
  TDirectory * newdir=NULL;
  TDirectory * chk = (TDirectory *) dir->Get("RunAverages"); 
  if(chk==NULL)  {
    newdir = (TDirectory *) dir->mkdir("RunAverages");
  } else {
    newdir = chk;
  }
  newdir->cd();
  r[roi].runcnt->Write(Form("cnt_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  r[roi].runqx->Write(Form("qx_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  r[roi].runqy->Write(Form("qy_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  savedir->cd();
}


bool Framework::AddFile(){
  char buf[120];
  if(fgets(buf,120,flist)==NULL) return false;
  buf[strlen(buf)-1]=0;
  string infile = buf;
  cout<<"infile:"<<infile<<":"<<endl;
  tf->Close();
  tf = new TFile(infile.data(),"read");
  tr = (TTree *) tf->Get("vnanalyzer/tree");
  maxevents = tr->GetEntries();
  tr->SetBranchAddress("Cent",&centval);
  tr->SetBranchAddress("NtrkOff",&Noff);
  tr->SetBranchAddress("ntrkflat",&ntrkval);
  tr->SetBranchAddress("Vtx",&vtx);
  tr->SetBranchAddress("epang",&epang);
  tr->SetBranchAddress("eporig",&eporig);
  tr->SetBranchAddress("qx",      &qx);
  tr->SetBranchAddress("qy",      &qy);
  tr->SetBranchAddress("q",       &q);
  tr->SetBranchAddress("vn", &vn );
  tr->SetBranchAddress("mult",    &epmult);
  tr->SetBranchAddress("sumw",    &sumw);
  tr->SetBranchAddress("sumw2",    &sumw2);
  tr->SetBranchAddress("Run",     &runno_);
  tr->SetBranchAddress("Rescor",  &rescor);
  tr->SetBranchAddress("RescorErr",&rescorErr);
  for(int iorder = 0; iorder<nqxorder; iorder++) {
    if(tr->FindBranch(Form("qxtrk%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qxtrk%d",qxorders[iorder]),  &qxtrk[iorder]);
    if(tr->FindBranch(Form("qytrk%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qytrk%d",qxorders[iorder]),  &qytrk[iorder]);
  }
  tr->SetBranchAddress("qcnt",    &qcnt);
  tr->SetBranchAddress("avpt",    &avpt);
  return true;
}
Framework::Framework(string filelist){
  ran = new TRandom();
  cout<<"open: "<<filelist<<endl;
  //system("cat filelist.dat");
  flist = fopen(filelist.data(),"r");
  char buf[120];
  fgets(buf,120,flist);
  buf[strlen(buf)-1]=0;
  string infile = buf;
  cout<<"infile:"<<infile<<":"<<endl;
  tf = new TFile(infile.data(),"read");
  tr = (TTree *) tf->Get("vnanalyzer/tree");
  maxevents = tr->GetEntries();
  tr->SetBranchAddress("Cent",&centval);
  tr->SetBranchAddress("NtrkOff",&Noff);
  tr->SetBranchAddress("ntrkflat",&ntrkval);
  tr->SetBranchAddress("Vtx",&vtx);
  tr->SetBranchAddress("epang",&epang);
  tr->SetBranchAddress("eporig",&eporig);
  tr->SetBranchAddress("qx",      &qx);
  tr->SetBranchAddress("qy",      &qy);
  tr->SetBranchAddress("q",       &q);
  tr->SetBranchAddress("vn", &vn );
  tr->SetBranchAddress("mult",    &epmult);
  tr->SetBranchAddress("sumw",    &sumw);
  tr->SetBranchAddress("sumw2",    &sumw2);
  tr->SetBranchAddress("Run",     &runno_);
  tr->SetBranchAddress("Rescor",  &rescor);
  tr->SetBranchAddress("RescorErr",&rescorErr);
  for(int iorder = 0; iorder<nqxorder; iorder++) {
    if(tr->FindBranch(Form("qxtrk%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qxtrk%d",qxorders[iorder]),  &qxtrk[iorder]);
    if(tr->FindBranch(Form("qytrk%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qytrk%d",qxorders[iorder]),  &qytrk[iorder]);
  }
  tr->SetBranchAddress("qcnt",    &qcnt);
  tr->SetBranchAddress("avpt",    &avpt);
  runs = new TH1D("runs","runs",maxRunRange-minRunRange+1,minRunRange,maxRunRange);
  runs->SetDirectory(0);
}
int Framework::SetROIRange(int order, int minCent, int maxCent, double minEta, double maxEta, double minPt, double maxPt) {
  if(maxevents==0) return -1;
  int ifound = -1;
  for(int i = 0; i<nqxorder; i++) {
    if(order == qxorders[i]) ifound = i;
  }
  if(ifound<0) {
    cout<<"The requested vn order is not part of the replay. ABORT!"<<endl;
    return -1;
  }
  r[nrange].orderIndx = ifound;
  GetEntry(0);
  r[nrange].order = order;
  r[nrange].minCent = minCent;
  r[nrange].maxCent = maxCent;
  r[nrange].minEtaBin = qxtrk[0]->GetYaxis()->FindBin(minEta+0.01);
  r[nrange].maxEtaBin = qxtrk[0]->GetYaxis()->FindBin(maxEta-0.01);
  r[nrange].minPtBin = qxtrk[0]->GetXaxis()->FindBin(minPt+0.01);
  r[nrange].maxPtBin = qxtrk[0]->GetXaxis()->FindBin(maxPt-0.01);
  r[nrange].minEta = qxtrk[0]->GetYaxis()->GetBinLowEdge(r[nrange].minEtaBin);
  r[nrange].maxEta = qxtrk[0]->GetYaxis()->GetBinLowEdge(r[nrange].maxEtaBin)+qxtrk[0]->GetYaxis()->GetBinWidth(r[nrange].maxEtaBin);
  r[nrange].minPt = qxtrk[0]->GetXaxis()->GetBinLowEdge(r[nrange].minPtBin);
  r[nrange].maxPt = qxtrk[0]->GetXaxis()->GetBinLowEdge(r[nrange].maxPtBin)+qxtrk[0]->GetXaxis()->GetBinWidth(r[nrange].maxPtBin);
  r[nrange].vn2d = new TH2D(Form("vn2d_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
			    Form("vn2d_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),100,-1.4,1.4,100,-1.4,1.4);
  r[nrange].vn2d->SetOption("colz");
  r[nrange].vn2d->SetDirectory(0);
  r[nrange].vn2d->SetXTitle(Form("v_{n,x}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn2d->SetYTitle("v_{n,y}^{obs}");
  r[nrange].vn1d = new TH1D(Form("vn1d_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
			    Form("vn1d_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),1000,0,1.4);
  r[nrange].vn1d->SetDirectory(0);
  r[nrange].vn1d->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn1d->SetYTitle("Counts");
  r[nrange].vn1dMult = new TH1D(Form("vn1dMult_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
			    Form("vn1dMult_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),1000,0,1.4);
  r[nrange].vn1dMult->SetDirectory(0);
  r[nrange].vn1dMult->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn1dMult->SetYTitle("AvMult");

  r[nrange].mult = new TH1D(Form("mult_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
			    Form("mult_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),120,0,120);
  r[nrange].mult->SetDirectory(0);
  if(runs!= NULL) {
    r[nrange].runcnt = new TH1D(Form("runcnt_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
				Form("runcnt_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runcnt->SetDirectory(0);
    r[nrange].runcnt->Sumw2();
    r[nrange].runqx = new TH1D(Form("runqx_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
			       Form("runqx_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runqx->SetDirectory(0);
    r[nrange].runqx->Sumw2();
    r[nrange].runqy = new TH1D(Form("runqy_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
			       Form("runqy_%d_%d_%d_%03.1f_%03.1f_%03.1f_%03.1f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runqy->SetDirectory(0);
    r[nrange].runqy->Sumw2();
  }
  if(!foff->IsZombie()) {
    r[nrange].qoffx = (TH1D *) foff->Get(Form("%d_%d/%3.1f_%3.1f/RunAverages/qx_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta));
    r[nrange].qoffx->Divide((TH1D *) foff->Get(Form("%d_%d/%3.1f_%3.1f/RunAverages/cnt_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta)));
    r[nrange].qoffx->SetDirectory(0);
    r[nrange].qoffy = (TH1D *)foff->Get(Form("%d_%d/%3.1f_%3.1f/RunAverages/qy_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta));
    r[nrange].qoffy->Divide((TH1D *) foff->Get(Form("%d_%d/%3.1f_%3.1f/RunAverages/cnt_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta)));
    r[nrange].qoffy->SetDirectory(0);

    r[nrange].qdifx = (TH1D *) foff->Get(Form("%d_%d/%3.1f_%3.1f/qdifx",minCent,maxCent,minPt,maxPt));
    r[nrange].qdify = (TH1D *) foff->Get(Form("%d_%d/%3.1f_%3.1f/qdify",minCent,maxCent,minPt,maxPt));
    r[nrange].qdifx->Divide((TH1D *) foff->Get(Form("%d_%d/%3.1f_%3.1f/qdifcnt",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdify->Divide((TH1D *) foff->Get(Form("%d_%d/%3.1f_%3.1f/qdifcnt",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdifx->SetDirectory(0);
    r[nrange].qdify->SetDirectory(0);

  }
  ++nrange;
  return nrange-1;
} 

void Framework::ShowAllROIRanges(){
  cout<<setw(12)<<right<<"indx"<<setw(12)<<right<<"minCent"<<setw(12)<<right<<"maxCent"<<setw(12)<<right<<"minEta"<<setw(12)<<right<<"maxEta"<<setw(12)<<right<<"minPt"<<setw(12)<<right<<"maxPt"<<setw(12)<<right<<"minEtaBin"<<setw(12)<<right<<"maxEtaBin"<<setw(12)<<right<<"minPtBin"<<setw(12)<<right<<"maxPtBin"<<std::endl;
  for(int i = 0; i<nrange; i++) {
    cout<<setprecision(2)<<setw(12)<<right<<i<<setprecision(2)<<setw(12)<<right<<r[i].minCent<<setprecision(2)<<setw(12)<<right<<r[i].maxCent<<setprecision(2)<<setw(12)<<right<<r[i].minEta<<setprecision(2)<<setw(12)<<right<<r[i].maxEta<<setprecision(2)<<setw(12)<<right<<r[i].minPt<<setprecision(2)<<setw(12)<<right<<r[i].maxPt<<setprecision(2)<<setw(12)<<right<<r[i].minEtaBin<<setprecision(2)<<setw(12)<<right<<r[i].maxEtaBin<<setprecision(2)<<setw(12)<<right<<r[i].minPtBin<<setprecision(2)<<setw(12)<<right<<r[i].maxPtBin<<std::endl;
  }
}

void Framework::SetROIEP(int roi, int EPA, int EPB, int EPC, int EPA2, int EPB2, int EPC2){
  if(roi>nrange) {
    cout<<"request roi index = "<<roi<< " > nrange = "<<nrange<<endl;
    return;
  }
  r[roi].A = EPA;
  r[roi].B = EPB;
  r[roi].C = EPC;
  r[roi].A2 = EPA2;
  r[roi].B2 = EPB2;
  r[roi].C2 = EPC2;
}

void Framework::AddEvent(int evt) {
  GetEntry(evt);
  runs->Fill((double) runno_);
  if(minRun>0 && runno_<minRun) return;
  if(maxRun>0 && runno_>maxRun) return;
  for(int i = 0; i<nrange; i++) {
    if(r[i].A<0) {
      cout<<"Event planes have not yet been initialized for ROI # "<<i<<endl;
    }
    if(centval<(double) r[i].minCent || centval>=(double)r[i].maxCent) continue;
    double qAx = qx[r[i].A];
    double qAy = qy[r[i].A];
    double qnx = 0;
    double qny = 0;
    double qnxe = 0;
    double qnye = 0;
    double qncnt = 0;
    double qncnte = 0;
    double pt = 0;
    for(int j = r[i].minEtaBin; j<= r[i].maxEtaBin; j++) {
      for(int k = r[i].minPtBin; k<=r[i].maxPtBin; k++) {
	qnx+=qxtrk[r[i].orderIndx]->GetBinContent(k,j);
	qny+=qytrk[r[i].orderIndx]->GetBinContent(k,j);
	qncnt+=qcnt->GetBinContent(k,j);
	if(qcnt->GetBinContent(k,j)>0) {
	  qnxe+=pow(qxtrk[r[i].orderIndx]->GetBinError(k,j),2);
	  qnye+=pow(qytrk[r[i].orderIndx]->GetBinError(k,j),2);
	  qncnte+=pow(qcnt->GetBinError(k,j),2);
	}
	pt+=avpt->GetBinContent(k,j);
      }
    }
    pt/=qncnt;
    r[i].mult->Fill(qncnt);
    vnxEvt[i] = -2;
    vnyEvt[i] = -2;
    if(minMult>0 && qncnt<minMult) continue;
    if(maxMult>0 && qncnt>maxMult) continue;
    int isub = ran->Uniform(0,9.9999);
    double val = qAx*qnx+qAy*qny;
    if(qncnt>0) {
      double qnxoff = 0;
      double qnyoff = 0;
      int runbin = 0;
      vnxEvt[i] = qnx/qncnt;
      vnyEvt[i] = qny/qncnt ;
      if(runs!=NULL) {
	runbin = runs->FindBin(runno_);
	if(!foff->IsZombie()){
	  qnxoff = r[i].qoffx->GetBinContent(runbin);
	  qnyoff = r[i].qoffy->GetBinContent(runbin);
	  vnxEvt[i]-=qnxoff;
	  vnyEvt[i]-=qnyoff;
	}
	r[i].runcnt->Fill(runno_);
	r[i].runqx->Fill(runno_,vnxEvt[i]);
	r[i].runqy->Fill(runno_,vnyEvt[i]);
      }
      r[i].vn2d->Fill(vnxEvt[i],vnyEvt[i]);
      r[i].vn1d->Fill(sqrt(pow(vnxEvt[i],2)+pow(vnyEvt[i],2)));
      r[i].vn1dMult->Fill(sqrt(pow(vnxEvt[i],2)+pow(vnyEvt[i],2)),qncnt);
      r[i].qn+=val;
      r[i].qne+=pow(qAx,2)*qnxe/fabs(qnx) + pow(qAy,2)*qnye/fabs(qny);
      r[i].qnSub[isub]+=val;
      r[i].qnSube[isub]+=pow(qAx,2)*qnxe/fabs(qnx) + pow(qAy,2)*qnye/fabs(qny);
      r[i].wn+=qncnt;
      r[i].wne+=qncnte;
      r[i].pt+=pt;
      double Ax = qx[r[i].A];
      double Bx = qx[r[i].B];
      double Cx = qx[r[i].C];
      double Ay = qy[r[i].A];
      double By = qy[r[i].B];
      double Cy = qy[r[i].C];
      double wA = sumw[r[i].A];
      double wB = sumw[r[i].B];
      double wC = sumw[r[i].C];
      r[i].qAB += Ax*Bx+Ay*By;
      r[i].qAC += Ax*Cx+Ay*Cy;
      r[i].qBC += Bx*Cx+By*Cy;
      r[i].wAB += wA*wB;
      r[i].wAC += wA*wC;
      r[i].wBC += wB*wC;
      r[i].wnA+=qncnt*wA;
      r[i].wnASub[isub]+=qncnt*wA;
    }
  }
}
double Framework::GetqABC(int roi) {
  double AB = r[roi].qAB/r[roi].wAB;
  double AC = r[roi].qAC/r[roi].wAC;
  double BC = r[roi].qBC/r[roi].wBC;
  return sqrt(AB*AC/BC);
}

TH1D * Framework::GetSpectra(int roi, bool effCorrect) {
  debug = false;
  TDirectory * dirsave = gDirectory;
  TH2D * spec;
  TH2D * heff=0;
  int minCent = r[roi].minCent;
  int maxCent = r[roi].maxCent;
  string crnge = Form("%d_%d",minCent,maxCent);
  string spname = "vnanalyzer/Spectra/"+crnge+"/ptspec";
  TH2D * ptcnt = (TH2D *) tf->Get(spname.data());
  spec = (TH2D *) ptcnt->Clone(Form("spec_%s",crnge.data()));
  spec->SetDirectory(0);
  double etamin = r[roi].minEta;
  double etamax = r[roi].maxEta;
  int ietamin = spec->GetYaxis()->FindBin(etamin+0.01);
  int ietamax = spec->GetYaxis()->FindBin(etamax-0.01);
  TH1D * spec1draw = (TH1D *) spec->ProjectionX(Form("spec1draw_%d_%d_%03.1f_%03.1f",minCent,maxCent,etamin,etamax),ietamin,ietamax);
  if(effCorrect){
    TFile * feff = new TFile(efffile.data());
    double avcent = (minCent+maxCent)/2.;
    if(avcent<5) {
      heff = (TH2D *) feff->Get("Eff_0_5")->Clone("Eff_clone");
    } else if (avcent<10) {
      heff = (TH2D *) feff->Get("Eff_5_10")->Clone("Eff_clone");
    } else if (avcent<30) {
      heff = (TH2D *) feff->Get("Eff_10_30")->Clone("Eff_clone");
    } else if (avcent<50) {
      heff = (TH2D *) feff->Get("Eff_30_50")->Clone("Eff_clone");
    } else if (avcent<100) {
      heff = (TH2D *) feff->Get("Eff_50_100")->Clone("Eff_clone");
    }
    heff->SetDirectory(0);
    for(int i = 1; i<=spec->GetNbinsX(); i++) {
      for(int j = 1; j<=spec->GetNbinsY(); j++) {
	double binpt = spec->GetXaxis()->GetBinCenter(i);
	double bineta = spec->GetYaxis()->GetBinCenter(j);
	double eff = heff->GetBinContent(heff->GetXaxis()->FindBin(bineta),heff->GetYaxis()->FindBin(binpt));
	if(eff == 0 && debug) cout<<i<<"\t"<<j<<"\t"<<eff<<"\t"<<binpt<<"\t"<<bineta<<endl;
	if(eff>0) spec->SetBinContent(i,j,spec->GetBinContent(i,j)/eff);
      }
    }
    feff->Close(); 
  }
  TH1D * spec1d = (TH1D *) spec->ProjectionX(Form("spec1d_%d_%d_%03.1f_%03.1f",minCent,maxCent,etamin,etamax),ietamin,ietamax);
  spec1d->SetDirectory(0);
  double bineta = 0.199;
  for(int i = 0; i<spec1d->GetNbinsX(); i++) {
    double deta = etamax-etamin;
    double dpt = spec1d->GetBinWidth(i+1);
    double binpt = spec1d->GetXaxis()->GetBinCenter(i+1);
    double eff = heff->GetBinContent(heff->GetXaxis()->FindBin(bineta),heff->GetYaxis()->FindBin(binpt));
    if(debug) cout<<spec1draw->GetXaxis()->GetBinCenter(i+1)<<"\t"<<spec1d->GetXaxis()->GetBinCenter(i+1)<<"\t"<<eff<<"\t"<<spec1draw->GetBinContent(i+1)<<"\t"<<spec1d->GetBinContent(i+1)<<"\t"<<deta<<"\t"<<dpt<<endl;
    spec1d->SetBinContent(i+1,spec1d->GetBinContent(i+1)/deta/dpt);
    spec1d->SetBinError(i+1,spec1d->GetBinError(i+1)/deta/dpt);

  }
  spec1d->SetXTitle("p_{T} (GeV/c)");
  spec1d->SetYTitle("dN/d#eta/dp_{T}");
  spec1d->SetMarkerStyle(20);
  spec1d->SetMarkerColor(kBlue);
  spec1d->SetLineColor(kBlue);
  spec->Delete();
  debug = false;
  dirsave->cd();
  return spec1d;
}

double Framework::GetVnErrSubEvt(int roi){
  double vnav=0;
  for(int i = 0; i<10; i++) vnav+=GetVnSub(roi,i);
  vnav/=10.;
  double sig = 0;
  for(int i = 0; i<10; i++) sig+=pow(GetVnSub(roi,i)-vnav,2);
  sig = sqrt(sig/9.);
  cout<<"err: "<<vnav<<"\t"<<sig<<endl;
  return sig;
}
#endif
