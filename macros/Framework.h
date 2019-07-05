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
  Framework(string filelist="filelist.dat", bool EffCorrect = true);
  bool AddFile();
  int GetN(){return maxevents;}
  void GetEntry(int i){tr->GetEntry(i);}
  int GetNrange(){return nrange;}
  int SetROIRange(int order, int minCent, int maxCent, double minEta, double maxEta, double minPt, double maxPt);
  void ShowAllROIRanges();
  void SetROIEP(int roi, int EPA, int EPB, int EPC, int EPA2=-1, int EPB2=-1, int EPC2 = -1);
  void AddEvent(int evt);
  double GetAvPt(int roi){return (r[roi].pt_p+r[roi].pt_m)/(r[roi].wn_p+r[roi].wn_m);}
  double GetAvPt_p(int roi){return r[roi].pt_p/r[roi].wn_p;}
  double GetAvPt_m(int roi){return r[roi].pt_m/r[roi].wn_m;}
  double GetqnA(int roi){return r[roi].qn/r[roi].wnA;}
  double GetqnA_p(int roi){return r[roi].qn_p/r[roi].wnA_p;}
  double GetqnA_m(int roi){return r[roi].qn_m/r[roi].wnA_m;}
  double GetVn(int roi) {return GetqnA(roi)/GetqABC(roi);}
  double GetVn_p(int roi) {return GetqnA_p(roi)/GetqABC(roi);}
  double GetVn_m(int roi) {return GetqnA_m(roi)/GetqABC(roi);}
  double GetAng(int roi) {return r[roi].ang;}
  double GetAng_p(int roi) {return r[roi].ang_p;}
  double GetAng_m(int roi) {return r[roi].ang_m;}
  double GetVnErrSubEvt(int roi);
  double GetVnErr(int roi) {return GetqnAError(roi)/GetqABC(roi);}
  double GetVnErr_p(int roi) {return GetqnAError_p(roi)/GetqABC(roi);}
  double GetVnErr_m(int roi) {return GetqnAError_m(roi)/GetqABC(roi);}
  double GetqABC(int roi);
  double GetAngHFp(int roi) {return r[roi].angHFp;}
  double GetAvEta(int roi) {return (r[roi].minEta+r[roi].maxEta)/2.;}
  TH1D * GetSpectra(int roi);
  TH1D * GetSpectra_p(int roi);
  TH1D * GetSpectra_m(int roi);
  TH2D * Get2d(int roi){return r[roi].vn2d;}
  TH2D * Get2d_p(int roi){return r[roi].vn2d_p;}
  TH2D * Get2d_m(int roi){return r[roi].vn2d_m;}
  TH1D * Get1d(int roi){return r[roi].vn1d;}
  TH1D * Get1d_p(int roi){return r[roi].vn1d_p;}
  TH1D * Get1d_m(int roi){return r[roi].vn1d_m;}
  TH1D * Get1dMult(int roi){return r[roi].vn1dMult;}
  TH1D * Get1dMult_p(int roi){return r[roi].vn1dMult_p;}
  TH1D * Get1dMult_m(int roi){return r[roi].vn1dMult_m;}
  int GetMinCent(int roi){return r[roi].minCent;}
  int GetMaxCent(int roi){return r[roi].maxCent;}
  int GetMinMult(){return minMult;}
  int GetMaxMult(){return maxMult;}
  int GetMinRun(){return minRun;}
  int GetMaxRun(){return maxRun;}
  bool Smear(int roi, double vnx,double vny, double &vnobsx, double &vnobsy, double &vnobs);
  bool Smear_p(int roi, double vnx,double vny, double &vnobsx, double &vnobsy, double &vnobs);
  bool Smear_m(int roi, double vnx,double vny, double &vnobsx, double &vnobsy, double &vnobs);
  TH1D * GetRuns(){return runs;}
  TH1D * GetMult(int roi){return r[roi].mult;}
  TH1D * GetMult_p(int roi){return r[roi].mult_p;}
  TH1D * GetMult_m(int roi){return r[roi].mult_m;}
  double GetVnxEvt(int roi){return vnxEvt[roi];}
  double GetVnyEvt(int roi){return vnyEvt[roi];}
  double GetVnxEvt_p(int roi){return vnxEvt_p[roi];}
  double GetVnyEvt_p(int roi){return vnyEvt_p[roi];}
  double GetVnxEvt_m(int roi){return vnxEvt_m[roi];}
  double GetVnyEvt_m(int roi){return vnyEvt_m[roi];}
  uint GetRunno(int evt){GetEntry(evt); return runno_;}
  void SetMinMult(int val){minMult = val;}
  void SetMaxMult(int val){maxMult = val;}
  void SetMinRun(int val){minRun = val;}
  void SetMaxRun(int val){maxRun = val;}
  void SetRuns(int nruns, double * runs); 
  void SaveFrameworkRuns(int roi, TDirectory * dir);
  void SaveFrameworkRuns_p(int roi, TDirectory * dir);
  void SaveFrameworkRuns_m(int roi, TDirectory * dir);
  bool LoadOffsets(string offname);
  bool LoadOffsets_p(string offname);
  bool LoadOffsets_m(string offname);
  double GetXdiff(int roi, int runno);
  double GetXdiff_p(int roi, int runno);
  double GetXdiff_m(int roi, int runno);
  double GetYdiff(int roi, int runno);
  double GetYdiff_p(int roi, int runno);
  double GetYdiff_m(int roi, int runno);
  string GetSide(int roi) {
    if(r[roi].minEta<0 && r[roi].maxEta<=0) {
      return "_neg_eta";
    } if( r[roi].minEta>=0 && r[roi].maxEta>0) {
      return "_pos_eta";
    } else {
      return "";
    }
  }
  FILE * flist;
private:
  TH1D * runbins=NULL;
  TRandom * ran;
  double vnxEvt[MAXROI];
  double vnyEvt[MAXROI];
  double vnxEvt_p[MAXROI];
  double vnyEvt_p[MAXROI];
  double vnxEvt_m[MAXROI];
  double vnyEvt_m[MAXROI];
  double GetVnSub(int roi,int i) {return r[roi].qnSub[i]/r[roi].wnASub[i]/GetqABC(roi);}
  double GetVnSub_p(int roi,int i) {return r[roi].qnSub_p[i]/r[roi].wnASub_p[i]/GetqABC(roi);}
  double GetVnSub_m(int roi,int i) {return r[roi].qnSub_m[i]/r[roi].wnASub_m[i]/GetqABC(roi);}
  double GetqnAError(int roi){return sqrt(r[roi].qne)/r[roi].wnA;}
  double GetqnAError_p(int roi){return sqrt(r[roi].qne_p)/r[roi].wnA_p;}
  double GetqnAError_m(int roi){return sqrt(r[roi].qne_m)/r[roi].wnA_m;}
  int maxevents;
  TFile * tf;
  TFile *foff;
  TFile *foff_p;
  TFile *foff_m;
  bool effCorrect = false;
  TFile * feff = NULL;
  TH2D * heff5 = 0;
  TH2D * heff10 = 0;
  TH2D * heff30 = 0;
  TH2D * heff50 = 0;
  TH2D * heff100 = 0;
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
  TH2F * qxtrk_p[7]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
  TH2F * qytrk_p[7]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
  TH2F * qxtrk_m[7]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
  TH2F * qytrk_m[7]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
  TH2F * qcnt=nullptr;
  TH2F * qcnt_p=nullptr;
  TH2F * qcnt_m=nullptr;
  TH2F * avpt=nullptr;
  TH2F * avpt_p=nullptr;
  TH2F * avpt_m=nullptr;
  TH1D * runs=NULL;
  TH1D * runcnt;
  TH1D * runcnt_p;
  TH1D * runcnt_m;
  TH1D * runqx;
  TH1D * runqx_p;
  TH1D * runqx_m;
  TH1D * runqy;
  TH1D * runqy_p;
  TH1D * runqy_m;
  TH1D * qdifx;
  TH1D * qdify;
  TH1D * qdifx_p;
  TH1D * qdify_p;
  TH1D * qdifx_m;
  TH1D * qdify_m;
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
    double qn_p=0;
    double qn_m=0;
    double qne=0;
    double qne_p=0;
    double qne_m=0;
    double pt = 0;
    double pt_p=0;
    double pt_m=0;
    double wn = 0;
    double wn_p=0;
    double wn_m=0;
    double wne=0;
    double wne_p=0;
    double wne_m=0;
    double wnA = 0;
    double wnA_p = 0;
    double wnA_m = 0;
    double qAB = 0;
    double wAB = 0;
    double qAC = 0;
    double wAC = 0;
    double qBC = 0;
    double wBC = 0;
    double qnSub[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSub_p[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSub_m[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSube[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSube_p[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSube_m[10]={0,0,0,0,0,0,0,0,0,0};
    double wnASub[10]={0,0,0,0,0,0,0,0,0,0};
    double wnASub_p[10]={0,0,0,0,0,0,0,0,0,0};
    double wnASub_m[10]={0,0,0,0,0,0,0,0,0,0};
    TH2D * vn2d;
    TH2D * vn2d_p;
    TH2D * vn2d_m;
    TH2D * smear;
    TH2D * smear_p;
    TH2D * smear_m;
    TH1D * vn1d;
    TH1D * vn1d_p;
    TH1D * vn1d_m;
    TH1D * vn1dMult;
    TH1D * vn1dMult_p;
    TH1D * vn1dMult_m;
    TH1D * mult;
    TH1D * mult_p;
    TH1D * mult_m;
    TH1D * runcnt;
    TH1D * runcnt_p;
    TH1D * runcnt_m;
    TH1D * runqx;
    TH1D * runqx_p;
    TH1D * runqx_m;
    TH1D * runqy;
    TH1D * runqy_p;
    TH1D * runqy_m;
    TH1D * qoffx;
    TH1D * qoffx_p;
    TH1D * qoffx_m;
    TH1D * qoffy;
    TH1D * qoffy_p;
    TH1D * qoffy_m;
    TH1D * qdifx;
    TH1D * qdifx_p;
    TH1D * qdifx_m;
    TH1D * qdify;
    TH1D * qdify_p;
    TH1D * qdify_m;
    double ang;
    double ang_p;
    double ang_m;
    double angHFp;
  } r[500];
  int nrange = 0;
};
bool Framework::LoadOffsets(string offname){
  foff = new TFile(offname.data(),"read"); 
  if(foff->IsZombie()) return false;

  return true;
}
bool Framework::LoadOffsets_p(string offname){
  foff_p = new TFile(offname.data(),"read"); 
  if(foff_p->IsZombie()) return false;

  return true;
}
bool Framework::LoadOffsets_m(string offname){
  foff_m = new TFile(offname.data(),"read"); 
  if(foff_m->IsZombie()) return false;

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
double Framework::GetXdiff_p(int roi, int runno){
  if(foff_p->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdifx_p->GetBinContent(runbin);
}
double Framework::GetXdiff_m(int roi, int runno){
  if(foff_m->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdifx_m->GetBinContent(runbin);
}
double Framework::GetYdiff(int roi, int runno){
  if(foff->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdify->GetBinContent(runbin);
}
double Framework::GetYdiff_p(int roi, int runno){
  if(foff_p->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdify_p->GetBinContent(runbin);
}
double Framework::GetYdiff_m(int roi, int runno){
  if(foff_m->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdify_m->GetBinContent(runbin);
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
void Framework::SaveFrameworkRuns_p(int roi, TDirectory * dir){
  TDirectory * savedir = gDirectory;
  TDirectory * newdir=NULL;
  TDirectory * chk = (TDirectory *) dir->Get("RunAverages"); 
  if(chk==NULL)  {
    newdir = (TDirectory *) dir->mkdir("RunAverages");
  } else {
    newdir = chk;
  }
  newdir->cd();
  r[roi].runcnt_p->Write(Form("cnt_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  r[roi].runqx_p->Write(Form("qx_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  r[roi].runqy_p->Write(Form("qy_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  savedir->cd();
}
void Framework::SaveFrameworkRuns_m(int roi, TDirectory * dir){
  TDirectory * savedir = gDirectory;
  TDirectory * newdir=NULL;
  TDirectory * chk = (TDirectory *) dir->Get("RunAverages"); 
  if(chk==NULL)  {
    newdir = (TDirectory *) dir->mkdir("RunAverages");
  } else {
    newdir = chk;
  }
  newdir->cd();
  r[roi].runcnt_m->Write(Form("cnt_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  r[roi].runqx_m->Write(Form("qx_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  r[roi].runqy_m->Write(Form("qy_%03.1f_%03.1f",r[roi].minEta,r[roi].maxEta));
  savedir->cd();
}

bool Framework::Smear(int roi, double vnx,double vny, double &vnobsx, double &vnobsy, double &vnobs){
  if(foff->IsZombie()) return false;
  double xs;
  double ys;
  r[roi].smear->GetRandom2(xs,ys);
  vnobsx = vnx+xs;
  vnobsy = vny+ys;
  vnobs = sqrt(pow(vnobsx,2)+pow(vnobsy,2));
  return true;
}
bool Framework::Smear_p(int roi, double vnx,double vny, double &vnobsx, double &vnobsy, double &vnobs){
  if(foff_p->IsZombie()) return false;
  double xs;
  double ys;
  r[roi].smear_p->GetRandom2(xs,ys);
  vnobsx = vnx+xs;
  vnobsy = vny+ys;
  vnobs = sqrt(pow(vnobsx,2)+pow(vnobsy,2));
  return true;
}
bool Framework::Smear_m(int roi, double vnx,double vny, double &vnobsx, double &vnobsy, double &vnobs){
  if(foff_m->IsZombie()) return false;
  double xs;
  double ys;
  r[roi].smear_m->GetRandom2(xs,ys);
  vnobsx = vnx+xs;
  vnobsy = vny+ys;
  vnobs = sqrt(pow(vnobsx,2)+pow(vnobsy,2));
  return true;
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
    if(tr->FindBranch(Form("qxtrk_p%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qxtrk_p%d",qxorders[iorder]),  &qxtrk_p[iorder]);
    if(tr->FindBranch(Form("qytrk_p%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qytrk_p%d",qxorders[iorder]),  &qytrk_p[iorder]);
    if(tr->FindBranch(Form("qxtrk_m%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qxtrk_m%d",qxorders[iorder]),  &qxtrk_m[iorder]);
    if(tr->FindBranch(Form("qytrk_m%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qytrk_m%d",qxorders[iorder]),  &qytrk_m[iorder]);
    if(qxtrk[iorder]==NULL) qxtrk[iorder] = (TH2F *) qxtrk_p[iorder]->Clone(Form("qxtrk_%d",iorder));
    if(qytrk[iorder]==NULL) qytrk[iorder] = (TH2F *) qytrk_p[iorder]->Clone(Form("qytrk_%d",iorder));
  }
  tr->SetBranchAddress("qcnt_p",    &qcnt_p);
  tr->SetBranchAddress("avpt_p",    &avpt_p);
  tr->SetBranchAddress("qcnt_m",    &qcnt_m);
  tr->SetBranchAddress("avpt_m",    &avpt_m);
  return true;
}
Framework::Framework(string filelist, bool EffCorrect){
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
    if(tr->FindBranch(Form("qxtrk_p%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qxtrk_p%d",qxorders[iorder]),  &qxtrk_p[iorder]);
    if(tr->FindBranch(Form("qytrk_p%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qytrk_p%d",qxorders[iorder]),  &qytrk_p[iorder]);
    if(tr->FindBranch(Form("qxtrk_m%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qxtrk_m%d",qxorders[iorder]),  &qxtrk_m[iorder]);
    if(tr->FindBranch(Form("qytrk_m%d",qxorders[iorder]))!=NULL) tr->SetBranchAddress(Form("qytrk_m%d",qxorders[iorder]),  &qytrk_m[iorder]);
    qxtrk[iorder]=NULL;
    qytrk[iorder]=NULL;
  }
  qcnt = NULL;
  avpt = NULL;
  tr->SetBranchAddress("qcnt_p",    &qcnt_p);
  tr->SetBranchAddress("avpt_p",    &avpt_p);
  tr->SetBranchAddress("qcnt_m",    &qcnt_m);
  tr->SetBranchAddress("avpt_m",    &avpt_m);
  runs = new TH1D("runs","runs",maxRunRange-minRunRange+1,minRunRange,maxRunRange);
  runs->SetDirectory(0);
  effCorrect = EffCorrect;
  if(effCorrect){
    feff = new TFile(efffile.data());
    heff5 = (TH2D *) feff->Get("Eff_0_5")->Clone("Eff_clone");
    heff10 = (TH2D *) feff->Get("Eff_5_10")->Clone("Eff_clone");
    heff30 = (TH2D *) feff->Get("Eff_10_30")->Clone("Eff_clone");
    heff50 = (TH2D *) feff->Get("Eff_30_50")->Clone("Eff_clone");
    heff100 = (TH2D *) feff->Get("Eff_50_100")->Clone("Eff_clone");
    heff5->SetDirectory(0);
    heff10->SetDirectory(0);
    heff30->SetDirectory(0);
    heff50->SetDirectory(0);
    heff100->SetDirectory(0);  
    feff->Close();  
  } else {
    cout<<"Efficienty correction is not being applied."<<endl;
  }
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
  r[nrange].minEtaBin = qxtrk_p[0]->GetYaxis()->FindBin(minEta+0.01);
  r[nrange].maxEtaBin = qxtrk_p[0]->GetYaxis()->FindBin(maxEta-0.01);
  r[nrange].minPtBin = qxtrk_p[0]->GetXaxis()->FindBin(minPt+0.01);
  r[nrange].maxPtBin = qxtrk_p[0]->GetXaxis()->FindBin(maxPt-0.01);
  r[nrange].minEta = qxtrk_p[0]->GetYaxis()->GetBinLowEdge(r[nrange].minEtaBin);
  r[nrange].maxEta = qxtrk_p[0]->GetYaxis()->GetBinLowEdge(r[nrange].maxEtaBin)+qxtrk_p[0]->GetYaxis()->GetBinWidth(r[nrange].maxEtaBin);
  r[nrange].minPt = qxtrk_p[0]->GetXaxis()->GetBinLowEdge(r[nrange].minPtBin);
  r[nrange].maxPt = qxtrk_p[0]->GetXaxis()->GetBinLowEdge(r[nrange].maxPtBin)+qxtrk_p[0]->GetXaxis()->GetBinWidth(r[nrange].maxPtBin);
  //
  r[nrange].vn2d = new TH2D(Form("vn2d_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("vn2d_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),100,-1.4,1.4,100,-1.4,1.4);
  r[nrange].vn2d->SetOption("colz");
  r[nrange].vn2d->SetDirectory(0);
  r[nrange].vn2d->SetXTitle(Form("v_{n,x}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn2d->SetYTitle("v_{n,y}^{obs}");

  r[nrange].vn2d_p = new TH2D(Form("vn2d_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("vn2d_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),100,-1.4,1.4,100,-1.4,1.4);
  r[nrange].vn2d_p->SetOption("colz");
  r[nrange].vn2d_p->SetDirectory(0);
  r[nrange].vn2d_p->SetXTitle(Form("v_{n,x}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn2d_p->SetYTitle("v_{n,y}^{obs}");

  r[nrange].vn2d_m = new TH2D(Form("vn2d_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("vn2d_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),100,-1.4,1.4,100,-1.4,1.4);
  r[nrange].vn2d_m->SetOption("colz");
  r[nrange].vn2d_m->SetDirectory(0);
  r[nrange].vn2d_m->SetXTitle(Form("v_{n,x}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn2d_m->SetYTitle("v_{n,y}^{obs}");

  //
  r[nrange].vn1d = new TH1D(Form("vn1d_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("vn1d_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),1000,0,1.4);
  r[nrange].vn1d->SetDirectory(0);
  r[nrange].vn1d->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn1d->SetYTitle("Counts");

  r[nrange].vn1d_p = new TH1D(Form("vn1d_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("vn1d_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),1000,0,1.4);
  r[nrange].vn1d_p->SetDirectory(0);
  r[nrange].vn1d_p->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn1d_p->SetYTitle("Counts");

  r[nrange].vn1d_m = new TH1D(Form("vn1d_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("vn1d_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),1000,0,1.4);
  r[nrange].vn1d_m->SetDirectory(0);
  r[nrange].vn1d_m->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn1d_m->SetYTitle("Counts");

  //
  r[nrange].vn1dMult = new TH1D(Form("vn1dMult_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  				Form("vn1dMult_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),1000,0,1.4);
  r[nrange].vn1dMult->SetDirectory(0);
  r[nrange].vn1dMult->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn1dMult->SetYTitle("AvMult");

  r[nrange].vn1dMult_p = new TH1D(Form("vn1dMult_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  				Form("vn1dMult_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),1000,0,1.4);
  r[nrange].vn1dMult_p->SetDirectory(0);
  r[nrange].vn1dMult_p->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn1dMult_p->SetYTitle("AvMult");

  r[nrange].vn1dMult_m = new TH1D(Form("vn1dMult_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  				Form("vn1dMult_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),1000,0,1.4);
  r[nrange].vn1dMult_m->SetDirectory(0);
  r[nrange].vn1dMult_m->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",minEta,maxEta));
  r[nrange].vn1dMult_m->SetYTitle("AvMult");

  //
  r[nrange].mult = new TH1D(Form("mult_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("mult_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),2000,0,2000);
  r[nrange].mult->SetDirectory(0);

  r[nrange].mult_p = new TH1D(Form("mult_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("mult_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),2000,0,2000);
  r[nrange].mult_p->SetDirectory(0);

  r[nrange].mult_m = new TH1D(Form("mult_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			    Form("mult_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),2000,0,2000);
  r[nrange].mult_m->SetDirectory(0);

  //
  if(runs!= NULL) {
    r[nrange].runcnt = new TH1D(Form("runcnt_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  				Form("runcnt_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runcnt->SetDirectory(0);
    r[nrange].runcnt->Sumw2();
    r[nrange].runqx = new TH1D(Form("runqx_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			       Form("runqx_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runqx->SetDirectory(0);
    r[nrange].runqx->Sumw2();
    r[nrange].runqy = new TH1D(Form("runqy_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			       Form("runqy_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runqy->SetDirectory(0);
    r[nrange].runqy->Sumw2();

    r[nrange].runcnt_p = new TH1D(Form("runcnt_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  				Form("runcnt_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runcnt_p->SetDirectory(0);
    r[nrange].runcnt_p->Sumw2();
    r[nrange].runqx_p = new TH1D(Form("runqx_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			       Form("runqx_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runqx_p->SetDirectory(0);
    r[nrange].runqx_p->Sumw2();
    r[nrange].runqy_p = new TH1D(Form("runqy_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			       Form("runqy_p_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runqy_p->SetDirectory(0);
    r[nrange].runqy_p->Sumw2();

    r[nrange].runcnt_m = new TH1D(Form("runcnt_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  				Form("runcnt_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runcnt_m->SetDirectory(0);
    r[nrange].runcnt_m->Sumw2();
    r[nrange].runqx_m = new TH1D(Form("runqx_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			       Form("runqx_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runqx_m->SetDirectory(0);
    r[nrange].runqx_m->Sumw2();
    r[nrange].runqy_m = new TH1D(Form("runqy_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),
  			       Form("runqy_m_%d_%d_%d_%03.1f_%03.1f_%04.2f_%04.2f",order,minCent,maxCent,minEta,maxEta,minPt,maxPt),nruns,runlist);
    r[nrange].runqy_m->SetDirectory(0);
    r[nrange].runqy_m->Sumw2();

  }
  /////////
  if(!foff->IsZombie()) {
    r[nrange].qoffx = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qx_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta));
    r[nrange].qoffx->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta)));
    r[nrange].qoffx->SetDirectory(0);

    r[nrange].qoffy = (TH1D *)foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qy_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta));
    r[nrange].qoffy->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta)));
    r[nrange].qoffy->SetDirectory(0);

    
    r[nrange].qdifx = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifx",minCent,maxCent,minPt,maxPt));
    r[nrange].qdify = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdify",minCent,maxCent,minPt,maxPt));
    r[nrange].qdifx->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdify->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdifx->SetDirectory(0);
    r[nrange].qdify->SetDirectory(0);    


    r[nrange].smear = (TH2D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/diff2d",minCent,maxCent,minPt,maxPt));
    r[nrange].smear->SetDirectory(0);
  }
  /////////
  if(!foff_p->IsZombie()) {
    r[nrange].qoffx_p = (TH1D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qx_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta));
    r[nrange].qoffx_p->Divide((TH1D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta)));
    r[nrange].qoffx_p->SetDirectory(0);

    r[nrange].qoffy_p = (TH1D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qy_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta));
    r[nrange].qoffy_p->Divide((TH1D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta)));
    r[nrange].qoffy_p->SetDirectory(0);

    r[nrange].qdifx_p = (TH1D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/qdifx",minCent,maxCent,minPt,maxPt));
    r[nrange].qdify_p = (TH1D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/qdify",minCent,maxCent,minPt,maxPt));
    r[nrange].qdifx_p->Divide((TH1D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdify_p->Divide((TH1D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdifx_p->SetDirectory(0);
    r[nrange].qdify_p->SetDirectory(0);    
    r[nrange].smear_p = (TH2D *) foff_p->Get(Form("%d_%d/%4.2f_%4.2f/diff2d",minCent,maxCent,minPt,maxPt));
    r[nrange].smear_p->SetDirectory(0);
  }
  /////////
  if(!foff_m->IsZombie()) {
    r[nrange].qoffx_m = (TH1D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qx_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta));
    r[nrange].qoffx_m->Divide((TH1D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta)));
    r[nrange].qoffx_m->SetDirectory(0);

    r[nrange].qoffy_m = (TH1D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qy_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta));
    r[nrange].qoffy_m->Divide((TH1D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_%3.1f_%3.1f",minCent,maxCent,minPt,maxPt,minEta,maxEta)));
    r[nrange].qoffy_m->SetDirectory(0);

    r[nrange].qdifx_m = (TH1D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/qdifx",minCent,maxCent,minPt,maxPt));
    r[nrange].qdify_m = (TH1D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/qdify",minCent,maxCent,minPt,maxPt));
    r[nrange].qdifx_m->Divide((TH1D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdify_m->Divide((TH1D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdifx_m->SetDirectory(0);
    r[nrange].qdify_m->SetDirectory(0);    
    r[nrange].smear_m = (TH2D *) foff_m->Get(Form("%d_%d/%4.2f_%4.2f/diff2d",minCent,maxCent,minPt,maxPt));
    r[nrange].smear_m->SetDirectory(0);
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
    r[i].angHFp = epang[HFp2];
    double qAx = qx[r[i].A];
    double qAy = qy[r[i].A];
    double qnx = 0;
    double qny = 0;
    double qncnt = 0;
    double qnxe = 0;
    double qnye = 0;
    double qncnte = 0;
    double qnx_p = 0;
    double qny_p = 0;
    double qnxe_p = 0;
    double qnye_p = 0;
    double qncnt_p = 0;
    double qncnte_p = 0;
    double qnx_m = 0;
    double qny_m = 0;
    double qnxe_m = 0;
    double qnye_m = 0;
    double qncnt_m = 0;
    double qncnte_m = 0;
    double pt = 0;
    double pt_p = 0;
    double pt_m = 0;
    if(qxtrk[r[i].orderIndx]==NULL) {
      qxtrk[r[i].orderIndx] = (TH2F *) qxtrk_p[r[i].orderIndx]->Clone(Form("qxtrk_%d",r[i].orderIndx));
      qxtrk[r[i].orderIndx]->Reset();
    }
    if(qytrk[r[i].orderIndx]==NULL) {
      qytrk[r[i].orderIndx] = (TH2F *) qytrk_p[r[i].orderIndx]->Clone(Form("qytrk_%d",r[i].orderIndx));
      qytrk[r[i].orderIndx]->Reset();
    }
    if(qcnt==NULL) {
      qcnt = (TH2F *) qcnt_p->Clone(Form("qcnt_%d",r[i].orderIndx));
      qcnt->Reset();
    }
    if(avpt==NULL) {
      avpt = (TH2F *) avpt_p->Clone(Form("avpt_%d",r[i].orderIndx));
      avpt->Reset();

    }
    for(int j = r[i].minEtaBin; j<= r[i].maxEtaBin; j++) {
      for(int k = r[i].minPtBin; k<=r[i].maxPtBin; k++) {
	double binpt = qxtrk_p[r[i].orderIndx]->GetXaxis()->GetBinCenter(k);
	double bineta = qxtrk_p[r[i].orderIndx]->GetYaxis()->GetBinCenter(j);
	double eff = 1;
	if(effCorrect) {
	  double avcent = centval;
	  TH2D * heff=0;
	  if(avcent<5) {
	    heff = heff5;
	  } else if (avcent<10) {
	    heff = heff10;
	  } else if (avcent<30) {
	    heff = heff30;
	  } else if (avcent<50) {
	    heff = heff50;
	  } else if (avcent<100) {
	    heff = heff100;
	  }
	  eff = heff->GetBinContent(heff->GetXaxis()->FindBin(bineta),heff->GetYaxis()->FindBin(binpt));
	}
	int ior = r[i].orderIndx;
	qxtrk[r[i].orderIndx]->Reset();
	qxtrk[r[i].orderIndx]->Add(qxtrk_p[r[i].orderIndx]);
	qxtrk[r[i].orderIndx]->Add(qxtrk_m[r[i].orderIndx]);

	qytrk[r[i].orderIndx]->Reset();
	qytrk[r[i].orderIndx]->Add(qytrk_p[r[i].orderIndx]);
	qytrk[r[i].orderIndx]->Add(qytrk_m[r[i].orderIndx]);

	qcnt->Reset();
	qcnt->Add(qcnt_p);
	qcnt->Add(qcnt_m);

	avpt->Reset();
	avpt->Add(avpt_p);
	avpt->Add(avpt_m);

	qnx+=qxtrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	qny+=qytrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	qncnt+=qcnt->GetBinContent(k,j)/eff;

	qnx_p+=qxtrk_p[r[i].orderIndx]->GetBinContent(k,j)/eff;
	qny_p+=qytrk_p[r[i].orderIndx]->GetBinContent(k,j)/eff;
	qncnt_p+=qcnt_p->GetBinContent(k,j)/eff;

	qnx_m+=qxtrk_m[r[i].orderIndx]->GetBinContent(k,j)/eff;
	qny_m+=qytrk_m[r[i].orderIndx]->GetBinContent(k,j)/eff;
	qncnt_m+=qcnt_m->GetBinContent(k,j)/eff;

	if(qcnt_p->GetBinContent(k,j)>0) {
	  qnxe_p+=pow(qxtrk_p[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qnye_p+=pow(qytrk_p[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qncnte_p+=pow(qcnt_p->GetBinError(k,j)/eff,2);
	}
	pt_p+=avpt_p->GetBinContent(k,j)/eff;

	if(qcnt_m->GetBinContent(k,j)>0) {
	  qnxe_m+=pow(qxtrk_m[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qnye_m+=pow(qytrk_m[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qncnte_m+=pow(qcnt_m->GetBinError(k,j)/eff,2);
	}
	pt_m+=avpt_m->GetBinContent(k,j)/eff;
	
	if(qcnt->GetBinContent(k,j)>0) {
	  qnxe+=pow(qxtrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qnye+=pow(qytrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qncnte+=pow(qcnt->GetBinError(k,j)/eff,2);
	}
	pt+=avpt->GetBinContent(k,j)/eff;
      }
    }
    
    r[i].mult->Fill(qncnt);
    r[i].mult_p->Fill(qncnt_p);
    r[i].mult_m->Fill(qncnt_m);
    int isub = ran->Uniform(0,9.9999);
    int runbin = 0;

    r[i].ang = -10;
    r[i].ang_p = -10;
    r[i].ang_m = -10;
    vnxEvt[i] = -10;
    vnxEvt_p[i] = -10;
    vnxEvt_m[i] = -10;
    vnyEvt[i] = -10;
    vnyEvt_p[i] = -10;
    vnyEvt_m[i] = -10;

    if(minMult<0 || qncnt>minMult) {
      if(maxMult<0 || qncnt<=maxMult) {
	double val = qAx*qnx+qAy*qny;
	double qnxoff = 0;
	double qnyoff = 0;
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
	r[i].ang = atan2(vnyEvt[i],vnxEvt[i])/r[i].order;
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
      }
    }
    
    if(minMult<0 || qncnt_p>minMult) {
      if(maxMult<0 || qncnt_p<=maxMult) {
	double val_p = qAx*qnx_p+qAy*qny_p;
	double qnxoff_p = 0;
	double qnyoff_p = 0;
	vnxEvt_p[i] = qnx_p/qncnt_p;
	vnyEvt_p[i] = qny_p/qncnt_p ;
	if(runs!=NULL) {
	  runbin = runs->FindBin(runno_);
	  if(!foff->IsZombie()){
	    qnxoff_p = r[i].qoffx_p->GetBinContent(runbin);
	    qnyoff_p = r[i].qoffy_p->GetBinContent(runbin);
	    vnxEvt_p[i]-=qnxoff_p;
	    vnyEvt_p[i]-=qnyoff_p;
	  }
	  r[i].runcnt_p->Fill(runno_);
	  r[i].runqx_p->Fill(runno_,vnxEvt_p[i]);
	  r[i].runqy_p->Fill(runno_,vnyEvt_p[i]);
	}
	r[i].ang_p = atan2(vnyEvt_p[i],vnxEvt_p[i])/r[i].order;
	r[i].vn2d_p->Fill(vnxEvt_p[i],vnyEvt_p[i]);
	r[i].vn1d_p->Fill(sqrt(pow(vnxEvt_p[i],2)+pow(vnyEvt_p[i],2)));
	r[i].vn1dMult_p->Fill(sqrt(pow(vnxEvt_p[i],2)+pow(vnyEvt_p[i],2)),qncnt_p);
	r[i].qn_p+=val_p;
	r[i].qne_p+=pow(qAx,2)*qnxe_p/fabs(qnx_p) + pow(qAy,2)*qnye_p/fabs(qny_p);
	r[i].qnSub_p[isub]+=val_p;
	r[i].qnSube_p[isub]+=pow(qAx,2)*qnxe_p/fabs(qnx_p) + pow(qAy,2)*qnye_p/fabs(qny_p);
	r[i].wn_p+=qncnt_p;
	r[i].wne_p+=qncnte_p;
	r[i].pt_p+=pt_p;
      }
    }

    if(minMult<0 || qncnt_m>minMult) {
      if(maxMult<0 || qncnt_m<=maxMult) {
	double val_m = qAx*qnx_m+qAy*qny_m;
	double qnxoff_m = 0;
	double qnyoff_m = 0;
	vnxEvt_m[i] = qnx_m/qncnt_m;
	vnyEvt_m[i] = qny_m/qncnt_m ;
	if(runs!=NULL) {
	  runbin = runs->FindBin(runno_);
	  if(!foff->IsZombie()){
	    qnxoff_m = r[i].qoffx_m->GetBinContent(runbin);
	    qnyoff_m = r[i].qoffy_m->GetBinContent(runbin);
	    vnxEvt_m[i]-=qnxoff_m;
	    vnyEvt_m[i]-=qnyoff_m;
	  }
	  r[i].runcnt_m->Fill(runno_);
	  r[i].runqx_m->Fill(runno_,vnxEvt_m[i]);
	  r[i].runqy_m->Fill(runno_,vnyEvt_m[i]);
	}    
	r[i].ang_m = atan2(vnyEvt_m[i],vnxEvt_m[i])/r[i].order;
	r[i].vn2d_m->Fill(vnxEvt_m[i],vnyEvt_m[i]);
	r[i].vn1d_m->Fill(sqrt(pow(vnxEvt_m[i],2)+pow(vnyEvt_m[i],2)));
	r[i].vn1dMult_m->Fill(sqrt(pow(vnxEvt_m[i],2)+pow(vnyEvt_m[i],2)),qncnt_m);
	r[i].qn_m+=val_m;
	r[i].qne_m+=pow(qAx,2)*qnxe_m/fabs(qnx_m) + pow(qAy,2)*qnye_m/fabs(qny_m);
	r[i].qnSub_m[isub]+=val_m;
	r[i].qnSube_m[isub]+=pow(qAx,2)*qnxe_m/fabs(qnx_m) + pow(qAy,2)*qnye_m/fabs(qny_m);
	r[i].wn_m+=qncnt_m;
	r[i].wne_m+=qncnte_m;
	r[i].pt_m+=pt_m;
      }
    }
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
    r[i].wnA_p+=qncnt_p*wA;
    r[i].wnA_m+=qncnt_m*wA;
    r[i].wnASub[isub]+=qncnt*wA;
    r[i].wnASub_p[isub]+=qncnt_p*wA;
    r[i].wnASub_m[isub]+=qncnt_m*wA;
  }
}

double Framework::GetqABC(int roi) {
  double AB = r[roi].qAB/r[roi].wAB;
  double AC = r[roi].qAC/r[roi].wAC;
  double BC = r[roi].qBC/r[roi].wBC;
  return sqrt(AB*AC/BC);
}

 TH1D * Framework::GetSpectra(int roi) {
  debug = false;
  TDirectory * dirsave = gDirectory;
  TH2D * spec;
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
  TH2D * heff=0;
  if(effCorrect){
    double avcent = (minCent+maxCent)/2.;
    if(avcent<5) {
      heff = heff5;
    } else if (avcent<10) {
      heff = heff10;
    } else if (avcent<30) {
      heff = heff30;
    } else if (avcent<50) {
      heff = heff50;
    } else if (avcent<100) {
      heff = heff100;
    }
    for(int i = 1; i<=spec->GetNbinsX(); i++) {
      for(int j = 1; j<=spec->GetNbinsY(); j++) {
	double binpt = spec->GetXaxis()->GetBinCenter(i);
	double bineta = spec->GetYaxis()->GetBinCenter(j);
	double eff = heff->GetBinContent(heff->GetXaxis()->FindBin(bineta),heff->GetYaxis()->FindBin(binpt));
	if(eff == 0 && debug) cout<<i<<"\t"<<j<<"\t"<<eff<<"\t"<<binpt<<"\t"<<bineta<<endl;
	if(eff>0) spec->SetBinContent(i,j,spec->GetBinContent(i,j)/eff);
      }
    } 
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
