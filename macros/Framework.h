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
void check(int i){cout<<"at location "<<i<<endl;}
static const int nqxorder = 1;
int qxorders[nqxorder]={2};
static const uint minRunRange = 326381;
static const uint maxRunRange = 327565;
static const int MAXROI = 500;
static const int MAXRUNS = 500; 
class Framework{
public:
  Framework(string filelist="filelist.dat", bool EffCorrect = true, int histbins = 100, double maxvn=0.6, int CS = 0);
  bool AddFile();
  int GetN(){return maxevents;}
  void GetEntry(int i){tr->GetEntry(i);}
  int GetNrange(){return nrange;}
  int SetROIRange(int order, int minCent, int maxCent, float *sub1, float *sub2, double minPt, double maxPt);
  void ShowAllROIRanges();
  void SetROIEP(int roi, int EPA1, int EPB1, int EPC1, int EPA2=-1, int EPB2=-1, int EPC2 = -1);  // separate for subevents 1 and 2
  void AddEvent(int evt);
  double GetAvPt(int roi){return (r[roi].pt/r[roi].wn_full);}
  double GetqnA_sub1(int roi){return r[roi].qn_sub1/r[roi].wnA_sub1;}
  double GetqnA_sub2(int roi){return r[roi].qn_sub2/r[roi].wnA_sub2;}
  double GetqnA_full(int roi){return r[roi].qn_full/r[roi].wnA_full;}
  double GetVn_sub1(int roi) {return GetqnA_sub1(roi)/GetqABC_sub1(roi);}
  double GetVn_sub2(int roi) {return GetqnA_sub2(roi)/GetqABC_sub2(roi);}
  double GetVn_full(int roi) {return GetqnA_full(roi)/GetqABC_full(roi);}
  double GetAng_sub1(int roi) {return r[roi].ang_sub1;}
  double GetAng_sub2(int roi) {return r[roi].ang_sub2;}
  double GetAng_full(int roi) {return r[roi].ang_full;}
  double GetVnErrSubEvt_sub1(int roi);
  double GetVnErrSubEvt_sub2(int roi);
  double GetVnErrSubEvt_full(int roi);
  double GetVnErr_sub1(int roi) {return GetqnAError_sub1(roi)/GetqABC_sub1(roi);}
  double GetVnErr_sub2(int roi) {return GetqnAError_sub2(roi)/GetqABC_sub2(roi);}
  double GetVnErr_full(int roi) {return GetqnAError_full(roi)/GetqABC_full(roi);}
  double GetqABC_sub1(int roi);
  double GetqABC_sub2(int roi);
  double GetqABC_full(int roi);
  double GetAngHFp(int roi) {return r[roi].angHFp;}
  double GetAvEta_sub1(int roi) {return (r[roi].minEta1+r[roi].maxEta1)/2.;}
  double GetAvEta_sub2(int roi) {return (r[roi].minEta2+r[roi].maxEta2)/2.;}
  double GetCent(int roi) {return r[roi].cent;}
  double GetVtx(int roi) {return r[roi].vtx;}
  TH1D * GetSpectra(int roi);
  TH2D * Get2d_sub1(int roi){return r[roi].vn2d_sub1;}
  TH2D * Get2d_sub2(int roi){return r[roi].vn2d_sub2;}
  TH2D * Get2d_full(int roi){return r[roi].vn2d_full;}
  TH1D * Get1d_sub1(int roi){return r[roi].vn1d_sub1;}
  TH1D * Get1d_sub2(int roi){return r[roi].vn1d_sub2;}
  TH1D * Get1d_full(int roi){return r[roi].vn1d_full;}
  TH1D * Get1dMult_sub1(int roi){return r[roi].vn1dMult_sub1;}
  TH1D * Get1dMult_sub2(int roi){return r[roi].vn1dMult_sub2;}
  TH1D * Get1dMult_full(int roi){return r[roi].vn1dMult_full;}
  int GetMinCent(int roi){return r[roi].minCent;}
  int GetMaxCent(int roi){return r[roi].maxCent;}
  int GetMinMult(){return minMult;}
  int GetMaxMult(){return maxMult;}
  int GetMinRun(){return minRun;}
  int GetMaxRun(){return maxRun;}
  bool Smear(int roi, double vnx,double vny, double &vnobsx, double &vnobsy, double &vnobs);
  TH1D * GetRuns(){return runs;}
  TH1D * GetMult_sub1(int roi){return r[roi].mult_sub1;}
  TH1D * GetMult_sub2(int roi){return r[roi].mult_sub2;}
  TH1D * GetMult_full(int roi){return r[roi].mult_full;}
  double GetVnxEvt_sub1(int roi){return vnxEvt_sub1[roi];}
  double GetVnyEvt_sub1(int roi){return vnyEvt_sub1[roi];}
  double GetVnxEvt_sub2(int roi){return vnxEvt_sub2[roi];}
  double GetVnyEvt_sub2(int roi){return vnyEvt_sub2[roi];}
  double GetVnxEvt_full(int roi){return vnxEvt_full[roi];}
  double GetVnyEvt_full(int roi){return vnyEvt_full[roi];}
  uint GetRunno(int evt){GetEntry(evt); return runno_;}
  void SetMinMult(int val){minMult = val;}
  void SetMaxMult(int val){maxMult = val;}
  void SetMinRun(int val){minRun = val;}
  void SetMaxRun(int val){maxRun = val;}
  void SetRuns(int nruns, double * runs); 
  void SaveFrameworkRuns(int roi, TDirectory * dir);
  bool LoadOffsets(string offname);
  double GetXdiff_sub1(int roi, int runno);
  double GetYdiff_sub1(int roi, int runno);
  double GetXdiff_sub2(int roi, int runno);
  double GetYdiff_sub2(int roi, int runno);
  double GetXdiff_full(int roi, int runno);
  double GetYdiff_full(int roi, int runno);
  FILE * flist;
private:
  int CS_;
  TH1D * runbins=NULL;
  TRandom * ran;
  double vnxEvt_sub1[MAXROI];
  double vnyEvt_sub1[MAXROI];
  double vnxEvt_sub2[MAXROI];
  double vnyEvt_sub2[MAXROI];
  double vnxEvt_full[MAXROI];
  double vnyEvt_full[MAXROI];
  double GetVnSub_sub1(int roi,int i) {return r[roi].qnSub_sub1[i]/r[roi].wnASub_sub1[i]/GetqABC_sub1(roi);}
  double GetVnSub_sub2(int roi,int i) {return r[roi].qnSub_sub2[i]/r[roi].wnASub_sub2[i]/GetqABC_sub2(roi);}
  double GetVnSub_full(int roi,int i) {return r[roi].qnSub_full[i]/r[roi].wnASub_full[i]/GetqABC_full(roi);}
  double GetqnAError_sub1(int roi){return sqrt(r[roi].qne_sub1)/r[roi].wnA_sub1;}
  double GetqnAError_sub2(int roi){return sqrt(r[roi].qne_sub2)/r[roi].wnA_sub2;}
  double GetqnAError_full(int roi){return sqrt(r[roi].qne_full)/r[roi].wnA_full;}
  int maxevents;
  TFile * tf;
  TFile *foff;
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
  TH1D * runqx;
  TH1D * runqy;
  TH1D * qdifx;
  TH1D * qdify;
  int nruns;
  double runlist[MAXRUNS];
  int nhistbins_;
  double maxvn_;
  struct range {
    int order;
    int orderIndx;
    int minCent;
    int maxCent;
    double minEta1;
    double maxEta1;
    double minEta2;
    double maxEta2;
    double minPt;
    double maxPt;
    int minEtaBin1;
    int maxEtaBin1;
    int minEtaBin2;
    int maxEtaBin2;
    int minPtBin;
    int maxPtBin;
    int A1=-1;
    int B1=-1;
    int C1=-1;
    int A2=-1;
    int B2=-1;
    int C2=-1;
    double qn_sub1=0;
    double qn_sub2=0;
    double qn_full=0;
    double qne_sub1=0;
    double qne_sub2=0;
    double qne_full=0;
    double pt = 0;
    double wn_sub1 = 0;
    double wn_sub2 = 0;
    double wn_full = 0;
    double wne_sub1=0;
    double wne_sub2=0;
    double wne_full=0;
    double wnA_sub1 = 0;
    double wnA_sub2 = 0;
    double wnA_full = 0;
    double qAB_sub1 = 0;
    double wAB_sub1 = 0;
    double qAC_sub1 = 0;
    double wAC_sub1 = 0;
    double qBC_sub1 = 0;
    double wBC_sub1 = 0;
    double qAB_sub2 = 0;
    double wAB_sub2 = 0;
    double qAC_sub2 = 0;
    double wAC_sub2 = 0;
    double qBC_sub2 = 0;
    double wBC_sub2 = 0;
    double qnSub_sub1[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSube_sub1[10]={0,0,0,0,0,0,0,0,0,0};
    double wnASub_sub1[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSub_sub2[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSube_sub2[10]={0,0,0,0,0,0,0,0,0,0};
    double wnASub_sub2[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSub_full[10]={0,0,0,0,0,0,0,0,0,0};
    double qnSube_full[10]={0,0,0,0,0,0,0,0,0,0};
    double wnASub_full[10]={0,0,0,0,0,0,0,0,0,0};
    TH2D * vn2d_sub1;
    TH2D * vn2d_sub2;
    TH2D * vn2d_full;
    TH2D * smear;
    TH1D * vn1d_sub1;
    TH1D * vn1d_sub2;
    TH1D * vn1d_full;
    TH1D * vn1dMult_sub1;
    TH1D * vn1dMult_sub2;
    TH1D * vn1dMult_full;
    TH1D * mult_sub1;
    TH1D * mult_sub2;
    TH1D * mult_full;
    TH1D * runcnt_sub1;
    TH1D * runqx_sub1;
    TH1D * runqy_sub1;
    TH1D * runcnt_sub2;
    TH1D * runqx_sub2;
    TH1D * runqy_sub2;
    TH1D * runcnt_full;
    TH1D * runqx_full;
    TH1D * runqy_full;
    TH1D * qoffx_sub1;
    TH1D * qoffy_sub1;
    TH1D * qoffx_sub2;
    TH1D * qoffy_sub2;
    TH1D * qoffx_full;
    TH1D * qoffy_full;
    TH1D * qdifx_sub1;
    TH1D * qdify_sub1;
    TH1D * qdifx_sub2;
    TH1D * qdify_sub2;
    TH1D * qdifx_full;
    TH1D * qdify_full;
    double ang_sub1;
    double ang_sub2;
    double ang_full;
    double vtx;
    double cent;
    double angHFp;
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

double Framework::GetXdiff_sub1(int roi, int runno){
  if(foff->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdifx_sub1->GetBinContent(runbin);
}
double Framework::GetYdiff_sub1(int roi, int runno){
  if(foff->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdify_sub1->GetBinContent(runbin);
}
double Framework::GetXdiff_sub2(int roi, int runno){
  if(foff->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdifx_sub2->GetBinContent(runbin);
}
double Framework::GetYdiff_sub2(int roi, int runno){
  if(foff->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdify_sub2->GetBinContent(runbin);
}
double Framework::GetXdiff_full(int roi, int runno){
  if(foff->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdifx_full->GetBinContent(runbin);
}
double Framework::GetYdiff_full(int roi, int runno){
  if(foff->IsZombie()) return 0;
  int runbin = runs->GetBin(runno);
  return r[roi].qdify_full->GetBinContent(runbin);
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
  r[roi].runcnt_sub1->Write("cnt_sub1");
  r[roi].runqx_sub1->Write("qx_sub1");
  r[roi].runqy_sub1->Write("qy_sub1");
  r[roi].runcnt_sub2->Write("cnt_sub2");
  r[roi].runqx_sub2->Write("qx_sub2");
  r[roi].runqy_sub2->Write("qy_sub2");
  r[roi].runcnt_full->Write("cnt_full");
  r[roi].runqx_full->Write("qx_full");
  r[roi].runqy_full->Write("qy_full");
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
Framework::Framework(string filelist, bool EffCorrect, int nhistbins, double maxvn, int CS){
  ran = new TRandom();
  nhistbins_ = nhistbins;
  maxvn_ = maxvn;
  CS_ = CS;
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
int Framework::SetROIRange(int order, int minCent, int maxCent, float *sub1, float *sub2, double minPt, double maxPt) {
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
  r[nrange].minEtaBin1 = qxtrk_p[0]->GetYaxis()->FindBin(sub1[0]+0.01);
  r[nrange].maxEtaBin1 = qxtrk_p[0]->GetYaxis()->FindBin(sub1[1]-0.01);
  r[nrange].minEtaBin2 = qxtrk_p[0]->GetYaxis()->FindBin(sub2[0]+0.01);
  r[nrange].maxEtaBin2 = qxtrk_p[0]->GetYaxis()->FindBin(sub2[1]-0.01);
  r[nrange].minPtBin = qxtrk_p[0]->GetXaxis()->FindBin(minPt+0.01);
  r[nrange].maxPtBin = qxtrk_p[0]->GetXaxis()->FindBin(maxPt-0.01);
  r[nrange].minEta1 = qxtrk_p[0]->GetYaxis()->GetBinLowEdge(r[nrange].minEtaBin1);
  r[nrange].maxEta1 = qxtrk_p[0]->GetYaxis()->GetBinLowEdge(r[nrange].maxEtaBin1)+qxtrk_p[0]->GetYaxis()->GetBinWidth(r[nrange].maxEtaBin1);
  r[nrange].minEta2 = qxtrk_p[0]->GetYaxis()->GetBinLowEdge(r[nrange].minEtaBin2);
  r[nrange].maxEta2 = qxtrk_p[0]->GetYaxis()->GetBinLowEdge(r[nrange].maxEtaBin2)+qxtrk_p[0]->GetYaxis()->GetBinWidth(r[nrange].maxEtaBin2);
  r[nrange].minPt = qxtrk_p[0]->GetXaxis()->GetBinLowEdge(r[nrange].minPtBin);
  r[nrange].maxPt = qxtrk_p[0]->GetXaxis()->GetBinLowEdge(r[nrange].maxPtBin)+qxtrk_p[0]->GetXaxis()->GetBinWidth(r[nrange].maxPtBin);
  //
  r[nrange].vn2d_sub1 = new TH2D(Form("vn2d_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  			         Form("vn2d_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),100,-1.4,1.4,100,-1.4,1.4);
  r[nrange].vn2d_sub1->SetOption("colz");
  r[nrange].vn2d_sub1->SetDirectory(0);
  r[nrange].vn2d_sub1->Sumw2();
  r[nrange].vn2d_sub1->SetXTitle(Form("v_{n,x}^{obs} (%3.1f < #eta <  %3.1f)",r[nrange].minEta1,r[nrange].maxEta1));
  r[nrange].vn2d_sub1->SetYTitle("v_{n,y}^{obs}");

  r[nrange].vn2d_sub2 = new TH2D(Form("vn2d_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),
  			         Form("vn2d_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),100,-1.4,1.4,100,-1.4,1.4);
  r[nrange].vn2d_sub2->SetOption("colz");
  r[nrange].vn2d_sub2->SetDirectory(0);
  r[nrange].vn2d_sub2->Sumw2();
  r[nrange].vn2d_sub2->SetXTitle(Form("v_{n,x}^{obs} (%3.1f < #eta <  %3.1f)",r[nrange].minEta2,r[nrange].maxEta2));
  r[nrange].vn2d_sub2->SetYTitle("v_{n,y}^{obs}");

  r[nrange].vn2d_full = new TH2D(Form("vn2d_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),
  			         Form("vn2d_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),100,-1.4,1.4,100,-1.4,1.4);
  r[nrange].vn2d_full->SetOption("colz");
  r[nrange].vn2d_full->SetDirectory(0);
  r[nrange].vn2d_full->Sumw2();
  r[nrange].vn2d_full->SetXTitle("v_{n,x}^{obs}");
  r[nrange].vn2d_full->SetYTitle("v_{n,y}^{obs}");


  //
  r[nrange].vn1d_sub1 = new TH1D(Form("vn1d_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  			         Form("vn1d_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),nhistbins_,0,maxvn_);
  r[nrange].vn1d_sub1->SetDirectory(0);
  r[nrange].vn1d_sub1->Sumw2();
  r[nrange].vn1d_sub1->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",r[nrange].minEta1,r[nrange].maxEta1));
  r[nrange].vn1d_sub1->SetYTitle("Counts");

  r[nrange].vn1d_sub2 = new TH1D(Form("vn1d_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),
  			         Form("vn1d_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),nhistbins_,0,maxvn_);
  r[nrange].vn1d_sub2->SetDirectory(0);
  r[nrange].vn1d_sub2->Sumw2();
  r[nrange].vn1d_sub2->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",r[nrange].minEta2,r[nrange].maxEta2));
  r[nrange].vn1d_sub2->SetYTitle("Counts");

  r[nrange].vn1d_full = new TH1D(Form("vn1d_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),
  			         Form("vn1d_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),nhistbins_,0,maxvn_);
  r[nrange].vn1d_full->SetDirectory(0);
  r[nrange].vn1d_full->Sumw2();
  r[nrange].vn1d_full->SetXTitle("v_{n}^{obs}");
  r[nrange].vn1d_full->SetYTitle("Counts");
  //
  r[nrange].vn1dMult_sub1 = new TH1D(Form("vn1dMult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  				     Form("vn1dMult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),nhistbins_,0,maxvn_);
  r[nrange].vn1dMult_sub1->SetDirectory(0);
  r[nrange].vn1dMult_sub1->Sumw2();
  r[nrange].vn1dMult_sub1->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",r[nrange].minEta1,r[nrange].maxEta1));
  r[nrange].vn1dMult_sub1->SetYTitle("AvMult");

  r[nrange].vn1dMult_sub2 = new TH1D(Form("vn1dMult_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),
  				     Form("vn1dMult_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),nhistbins_,0,maxvn_);
  r[nrange].vn1dMult_sub2->SetDirectory(0);
  r[nrange].vn1dMult_sub2->Sumw2();
  r[nrange].vn1dMult_sub2->SetXTitle(Form("v_{n}^{obs} (%3.1f < #eta <  %3.1f)",r[nrange].minEta2,r[nrange].maxEta2));
  r[nrange].vn1dMult_sub2->SetYTitle("AvMult");

  r[nrange].vn1dMult_full = new TH1D(Form("vn1dMult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  				     Form("vn1dMult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),nhistbins_,0,maxvn_);
  r[nrange].vn1dMult_full->SetDirectory(0);
  r[nrange].vn1dMult_full->Sumw2();
  r[nrange].vn1dMult_full->SetXTitle("v_{n}^{obs}");
  r[nrange].vn1dMult_full->SetYTitle("AvMult");


  //
  r[nrange].mult_sub1 = new TH1D(Form("mult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  			         Form("mult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),2000,0,2000);
  r[nrange].mult_sub1->SetDirectory(0);
  r[nrange].mult_sub1->Sumw2();

  r[nrange].mult_sub2 = new TH1D(Form("mult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  			         Form("mult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),2000,0,2000);
  r[nrange].mult_sub2->SetDirectory(0);
  r[nrange].mult_sub2->Sumw2();

  r[nrange].mult_full = new TH1D(Form("mult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  			         Form("mult_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),2000,0,2000);
  r[nrange].mult_full->SetDirectory(0);
  r[nrange].mult_full->Sumw2();
  //
  if(runs!= NULL) {
    r[nrange].runcnt_sub1 = new TH1D(Form("runcnt_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  				     Form("runcnt_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runcnt_sub1->SetDirectory(0);
    r[nrange].runcnt_sub1->Sumw2();
    r[nrange].runqx_sub1 = new TH1D(Form("runqx_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  			            Form("runqx_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runqx_sub1->SetDirectory(0);
    r[nrange].runqx_sub1->Sumw2();
    r[nrange].runqy_sub1 = new TH1D(Form("runqy_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),
  			            Form("runqy_%d_%d_%d_%04.2f_%04.2f_sub1",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runqy_sub1->SetDirectory(0);
    r[nrange].runqy_sub1->Sumw2();

    r[nrange].runcnt_sub2 = new TH1D(Form("runcnt_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),
  				     Form("runcnt_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runcnt_sub2->SetDirectory(0);
    r[nrange].runcnt_sub2->Sumw2();
    r[nrange].runqx_sub2 = new TH1D(Form("runqx_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),
  			            Form("runqx_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runqx_sub2->SetDirectory(0);
    r[nrange].runqx_sub2->Sumw2();
    r[nrange].runqy_sub2 = new TH1D(Form("runqy_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),
  			            Form("runqy_%d_%d_%d_%04.2f_%04.2f_sub2",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runqy_sub2->SetDirectory(0);
    r[nrange].runqy_sub2->Sumw2();

    r[nrange].runcnt_full = new TH1D(Form("runcnt_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),
  				     Form("runcnt_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runcnt_full->SetDirectory(0);
    r[nrange].runcnt_full->Sumw2();
    r[nrange].runqx_full = new TH1D(Form("runqx_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),
  			            Form("runqx_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runqx_full->SetDirectory(0);
    r[nrange].runqx_full->Sumw2();
    r[nrange].runqy_full = new TH1D(Form("runqy_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),
  			            Form("runqy_%d_%d_%d_%04.2f_%04.2f_full",order,minCent,maxCent,minPt,maxPt),nruns,runlist);
    r[nrange].runqy_full->SetDirectory(0);
    r[nrange].runqy_full->Sumw2();

  }
  /////////
  if(!foff->IsZombie()) {
    r[nrange].qoffx_sub1 = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qx_sub1",minCent,maxCent,minPt,maxPt));
    r[nrange].qoffx_sub1->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_sub1",minCent,maxCent,minPt,maxPt)));
    r[nrange].qoffx_sub1->SetDirectory(0);
    r[nrange].qoffy_sub1 = (TH1D *)foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qy_sub1",minCent,maxCent,minPt,maxPt));
    r[nrange].qoffy_sub1->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_sub1",minCent,maxCent,minPt,maxPt)));
    r[nrange].qoffy_sub1->SetDirectory(0);
    
    r[nrange].qoffx_sub2 = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qx_sub2",minCent,maxCent,minPt,maxPt));
    r[nrange].qoffx_sub2->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_sub2",minCent,maxCent,minPt,maxPt)));
    r[nrange].qoffx_sub2->SetDirectory(0);
    r[nrange].qoffy_sub2 = (TH1D *)foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qy_sub2",minCent,maxCent,minPt,maxPt));
    r[nrange].qoffy_sub2->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_sub2",minCent,maxCent,minPt,maxPt)));
    r[nrange].qoffy_sub2->SetDirectory(0);
    
    r[nrange].qoffx_full = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qx_full",minCent,maxCent,minPt,maxPt));
    r[nrange].qoffx_full->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_full",minCent,maxCent,minPt,maxPt)));
    r[nrange].qoffx_full->SetDirectory(0);
    r[nrange].qoffy_full = (TH1D *)foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/qy_full",minCent,maxCent,minPt,maxPt));
    r[nrange].qoffy_full->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/RunAverages/cnt_full",minCent,maxCent,minPt,maxPt)));
    r[nrange].qoffy_full->SetDirectory(0);
    
    r[nrange].qdifx_sub1 = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifx_sub1",minCent,maxCent,minPt,maxPt));
    r[nrange].qdify_sub1 = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdify_sub1",minCent,maxCent,minPt,maxPt));
    r[nrange].qdifx_sub1->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt_sub1",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdify_sub1->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt_sub1",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdifx_sub1->SetDirectory(0);
    r[nrange].qdify_sub1->SetDirectory(0);

    r[nrange].qdifx_sub2 = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifx_sub2",minCent,maxCent,minPt,maxPt));
    r[nrange].qdify_sub2 = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdify_sub12",minCent,maxCent,minPt,maxPt));
    r[nrange].qdifx_sub2->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt_sub2",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdify_sub2->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt_sub2",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdifx_sub2->SetDirectory(0);
    r[nrange].qdify_sub2->SetDirectory(0);

    r[nrange].qdifx_full = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifx_full",minCent,maxCent,minPt,maxPt));
    r[nrange].qdify_full = (TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdify_full",minCent,maxCent,minPt,maxPt));
    r[nrange].qdifx_full->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt_full",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdify_full->Divide((TH1D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/qdifcnt_full",minCent,maxCent,minPt,maxPt)));
    r[nrange].qdifx_full->SetDirectory(0);
    r[nrange].qdify_full->SetDirectory(0);


    r[nrange].smear = (TH2D *) foff->Get(Form("%d_%d/%4.2f_%4.2f/diff2d",minCent,maxCent,minPt,maxPt));
    r[nrange].smear->SetDirectory(0);
  }
  
  ++nrange;
  return nrange-1;
} 

void Framework::ShowAllROIRanges(){
  cout<<setw(12)<<right<<"indx"<<setw(12)<<right<<"minCent"<<setw(12)<<right<<"maxCent"<<setw(12)<<right<<"minEta1"<<setw(12)<<right<<"maxEta1"<<setw(12)<<right<<"minEta2"<<setw(12)<<right<<"maxEta2"<<setw(12)<<right<<"minPt"<<setw(12)<<right<<"maxPt"<<setw(12)<<right<<"minEtaBin1"<<setw(12)<<right<<"maxEtaBin1"<<setw(12)<<right<<"minEtaBin2"<<setw(12)<<right<<"maxEtaBin2"<<setw(12)<<right<<"minPtBin"<<setw(12)<<right<<"maxPtBin"<<std::endl;
  for(int i = 0; i<nrange; i++) {
    cout<<setprecision(2)<<setw(12)<<right<<i<<setprecision(2)<<setw(12)<<right<<r[i].minCent<<setprecision(2)<<setw(12)<<right<<r[i].maxCent<<setprecision(2)<<setw(12)<<right<<r[i].minEta1<<setprecision(2)<<setw(12)<<right<<r[i].maxEta1<<setw(12)<<right<<r[i].minEta2<<setprecision(2)<<setw(12)<<right<<r[i].maxEta2<<setprecision(2)<<setw(12)<<right<<r[i].minPt<<setprecision(2)<<setw(12)<<right<<r[i].maxPt<<setprecision(2)<<setw(12)<<right<<r[i].minEtaBin1<<setprecision(2)<<setw(12)<<right<<r[i].maxEtaBin1<<setw(12)<<right<<r[i].minEtaBin2<<setprecision(2)<<setw(12)<<right<<r[i].maxEtaBin2<<setprecision(2)<<setw(12)<<right<<r[i].minPtBin<<setprecision(2)<<setw(12)<<right<<r[i].maxPtBin<<std::endl;
  }
}

void Framework::SetROIEP(int roi, int EPA1, int EPB1, int EPC1, int EPA2, int EPB2, int EPC2){
  if(roi>nrange) {
    cout<<"request roi index = "<<roi<< " > nrange = "<<nrange<<endl;
    return;
  }
  r[roi].A1 = EPA1;
  r[roi].B1 = EPB1;
  r[roi].C1 = EPC1;
  r[roi].A2 = EPA2;
  r[roi].B2 = EPB2;
  r[roi].C2 = EPC2;
}

void Framework::AddEvent(int evt) {
  GetEntry(evt);
  runs->Fill((double) runno_);
  if(minRun>0 && runno_<minRun) return;
  if(maxRun>0 && runno_>maxRun) return;
  for(int iorder = 0; iorder<nqxorder; iorder++) {
    if(qxtrk[iorder]==NULL) {
      if(CS_ >=0){
	qxtrk[iorder] = (TH2F *) qxtrk_p[iorder]->Clone(Form("qxtrk_%d",iorder));
	if(CS_==0) qxtrk[iorder]->Add(qxtrk_m[iorder]);
      } else {
	qxtrk[iorder] = (TH2F *) qxtrk_m[iorder]->Clone(Form("qxtrk_%d",iorder));
      }
    } else {
      qxtrk[iorder]->Reset();
      if(CS_ >=0){
	qxtrk[iorder]->Add(qxtrk_p[iorder]);
	if(CS_==0) qxtrk[iorder]->Add(qxtrk_m[iorder]);
      } else {
	qxtrk[iorder]->Add(qxtrk_m[iorder]);
      }
    }
    if(qytrk[iorder]==NULL) {
      if(CS_ >=0) {
	qytrk[iorder] = (TH2F *) qytrk_p[iorder]->Clone(Form("qytrk_%d",iorder));
	if(CS_==0) {
	  qytrk[iorder]->Add(qytrk_m[iorder]);
	}
      } else {
	qytrk[iorder] = (TH2F *) qytrk_m[iorder]->Clone(Form("qytrk_%d",iorder));
      }
    } else{
      qytrk[iorder]->Reset();
      if(CS_ >=0){
	qytrk[iorder]->Add(qytrk_p[iorder]);
	if(CS_==0) qytrk[iorder]->Add(qytrk_m[iorder]);
      } else {
	qytrk[iorder]->Add(qytrk_m[iorder]);
      }
    }
  }
  if(qcnt==NULL) {
    if(CS_>=0) {
      qcnt = (TH2F *) qcnt_p->Clone("qcnt");
      if(CS_==0) qcnt->Add(qcnt_m);
    } else {
      qcnt = (TH2F *) qcnt_m->Clone("qcnt");
    }
  } else {
    qcnt->Reset();
    if(CS_>=0) {
      qcnt->Add( qcnt_p);
      if(CS_==0) qcnt->Add(qcnt_m);
    } else {
      qcnt->Add(qcnt_m);
    }
  }
  if(avpt==NULL) {
    if(CS_>=0) {
      avpt = (TH2F *) avpt_p->Clone("avpt"); 
      if(CS_==0) avpt->Add(avpt_m);
    } else {
      avpt = (TH2F *) avpt_m->Clone("avpt"); 
    }   
  } else {
    avpt->Reset();
    if(CS_>=0) {
      avpt->Add(avpt_p);
      if(CS_==0) avpt->Add(avpt_m);
    } else {
      avpt->Add(avpt_m);
    }
  }
  

  for(int i = 0; i<nrange; i++) {
    if(r[i].A1<0) {
      cout<<"Event planes have not yet been initialized for ROI # "<<i<<endl;
    }
    if(centval<(double) r[i].minCent || centval>=(double)r[i].maxCent) continue;
    r[i].angHFp = epang[HFp2];
    r[i].cent = centval;
    r[i].vtx = vtx;
    double qAx_sub1 = qx[r[i].A1];
    double qAy_sub1 = qy[r[i].A1];
    double qAx_sub2 = qx[r[i].A2];
    double qAy_sub2 = qy[r[i].A2];
    double qAx_full = qx[r[i].A1] + qx[r[i].A2];
    double qAy_full = qy[r[i].A2] + qy[r[i].A2];;
    double qnx_sub1 = 0;
    double qny_sub1 = 0;
    double qncnt_sub1 = 0;
    double qnxe_sub1 = 0;
    double qnye_sub1 = 0;
    double qncnte_sub1 = 0;

    double qnx_sub2 = 0;
    double qny_sub2 = 0;
    double qncnt_sub2 = 0;
    double qnxe_sub2 = 0;
    double qnye_sub2 = 0;
    double qncnte_sub2 = 0;

    double qnx_full = 0;
    double qny_full = 0;
    double qncnt_full = 0;
    double qnxe_full = 0;
    double qnye_full = 0;
    double qncnte_full = 0;
    double pt = 0;
    // sub1
    for(int j = r[i].minEtaBin1; j<= r[i].maxEtaBin1; j++) {
      for(int k = r[i].minPtBin; k<=r[i].maxPtBin; k++) {
	double binpt = qxtrk[r[i].orderIndx]->GetXaxis()->GetBinCenter(k);
	double bineta = qxtrk[r[i].orderIndx]->GetYaxis()->GetBinCenter(j);
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

	if(qcnt->GetBinContent(k,j)>0) {
	  qnx_sub1+=qxtrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	  qny_sub1+=qytrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	  qncnt_sub1+=qcnt->GetBinContent(k,j)/eff;
	
	  qnxe_sub1+=pow(qxtrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qnye_sub1+=pow(qytrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qncnte_sub1+=pow(qcnt->GetBinError(k,j)/eff,2);

	  qnx_full+=qxtrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	  qny_full+=qytrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	  qncnt_full+=qcnt->GetBinContent(k,j)/eff;
	
	  qnxe_full+=pow(qxtrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qnye_full+=pow(qytrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qncnte_full+=pow(qcnt->GetBinError(k,j)/eff,2);
	}
	pt+=avpt->GetBinContent(k,j)/eff;
      }
    }
    r[i].mult_sub1->Fill(qncnt_sub1);
    // sub2
    for(int j = r[i].minEtaBin2; j<= r[i].maxEtaBin2; j++) {
      for(int k = r[i].minPtBin; k<=r[i].maxPtBin; k++) {
	double binpt = qxtrk[r[i].orderIndx]->GetXaxis()->GetBinCenter(k);
	double bineta = qxtrk[r[i].orderIndx]->GetYaxis()->GetBinCenter(j);
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

	if(qcnt->GetBinContent(k,j)>0) {
	  qnx_sub2+=qxtrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	  qny_sub2+=qytrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	  qncnt_sub2+=qcnt->GetBinContent(k,j)/eff;
	
	  qnxe_sub2+=pow(qxtrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qnye_sub2+=pow(qytrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qncnte_sub2+=pow(qcnt->GetBinError(k,j)/eff,2);

	  qnx_full+=qxtrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	  qny_full+=qytrk[r[i].orderIndx]->GetBinContent(k,j)/eff;
	  qncnt_full+=qcnt->GetBinContent(k,j)/eff;
	
	  qnxe_full+=pow(qxtrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qnye_full+=pow(qytrk[r[i].orderIndx]->GetBinError(k,j)/eff,2);
	  qncnte_full+=pow(qcnt->GetBinError(k,j)/eff,2);
	}
	pt+=avpt->GetBinContent(k,j)/eff;
      }
    }
    r[i].mult_sub2->Fill(qncnt_sub2);
    r[i].mult_full->Fill(qncnt_full);

    int isub = ran->Uniform(0,9.9999);
    int runbin = 0;

    r[i].ang_sub1 = -10;
    r[i].ang_sub2 = -10;
    r[i].ang_full = -10;
    vnxEvt_sub1[i] = -10;
    vnyEvt_sub1[i] = -10;
    vnxEvt_sub2[i] = -10;
    vnyEvt_sub2[i] = -10;
    vnxEvt_full[i] = -10;
    vnyEvt_full[i] = -10;

    if(minMult<0 || (qncnt_sub1>minMult && qncnt_sub2>minMult)) {
      if(maxMult<0 || qncnt_full<=maxMult) {
	double val_sub1 = qAx_sub1*qnx_sub1+qAy_sub1*qny_sub1;
	double qnxoff_sub1 = 0;
	double qnyoff_sub1 = 0;
	double val_sub2 = qAx_sub2*qnx_sub2+qAy_sub2*qny_sub2;
	double qnxoff_sub2 = 0;
	double qnyoff_sub2 = 0;
	double val_full = qAx_full*qnx_full+qAy_full*qny_full;
	double qnxoff_full = 0;
	double qnyoff_full = 0;
	vnxEvt_sub1[i] = qnx_sub1/qncnt_sub1;
	vnxEvt_sub2[i] = qnx_sub2/qncnt_sub2;
	vnxEvt_full[i] = qnx_full/qncnt_full;
	vnyEvt_sub1[i] = qny_sub1/qncnt_sub1 ;
	vnyEvt_sub2[i] = qny_sub2/qncnt_sub2 ;
	vnyEvt_full[i] = qny_full/qncnt_full ;
	if(runs!=NULL) {
	  runbin = runs->FindBin(runno_);
	  if(!foff->IsZombie()){
	    qnxoff_sub1 = r[i].qoffx_sub1->GetBinContent(runbin);
	    qnyoff_sub1 = r[i].qoffy_sub1->GetBinContent(runbin);
	    vnxEvt_sub1[i]-=qnxoff_sub1;
	    vnyEvt_sub1[i]-=qnyoff_sub1;

	    qnxoff_sub2 = r[i].qoffx_sub2->GetBinContent(runbin);
	    qnyoff_sub2 = r[i].qoffy_sub2->GetBinContent(runbin);
	    vnxEvt_sub2[i]-=qnxoff_sub2;
	    vnyEvt_sub2[i]-=qnyoff_sub2;

	    qnxoff_full = r[i].qoffx_full->GetBinContent(runbin);
	    qnyoff_full = r[i].qoffy_full->GetBinContent(runbin);
	    vnxEvt_full[i]-=qnxoff_full;
	    vnyEvt_full[i]-=qnyoff_full;
	  }
	  r[i].runcnt_sub1->Fill(runno_);
	  r[i].runqx_sub1->Fill(runno_,vnxEvt_sub1[i]);
	  r[i].runqy_sub1->Fill(runno_,vnyEvt_sub1[i]);
	  r[i].runcnt_sub2->Fill(runno_);
	  r[i].runqx_sub2->Fill(runno_,vnxEvt_sub2[i]);
	  r[i].runqy_sub2->Fill(runno_,vnyEvt_sub2[i]);
	  r[i].runcnt_full->Fill(runno_);
	  r[i].runqx_full->Fill(runno_,vnxEvt_full[i]);
	  r[i].runqy_full->Fill(runno_,vnyEvt_full[i]);
	}
	r[i].ang_sub1 = atan2(vnyEvt_sub1[i],vnxEvt_sub1[i])/r[i].order;
	r[i].vn2d_sub1->Fill(vnxEvt_sub1[i],vnyEvt_sub1[i]);
	r[i].vn1d_sub1->Fill(sqrt(pow(vnxEvt_sub1[i],2)+pow(vnyEvt_sub1[i],2)));
	r[i].vn1dMult_sub1->Fill(sqrt(pow(vnxEvt_sub1[i],2)+pow(vnyEvt_sub1[i],2)),qncnt_sub1);
	r[i].qn_sub1+=val_sub1;
	r[i].qne_sub1+=pow(qAx_sub1,2)*qnxe_sub1/fabs(qnx_sub1) + pow(qAy_sub1,2)*qnye_sub1/fabs(qny_sub1);
	r[i].qnSub_sub1[isub]+=val_sub1;
	r[i].qnSube_sub1[isub]+=pow(qAx_sub1,2)*qnxe_sub1/fabs(qnx_sub1) + pow(qAy_sub1,2)*qnye_sub1/fabs(qny_sub1);
	r[i].wn_sub1+=qncnt_sub1;
	r[i].wne_sub1+=qncnte_sub1;

	r[i].ang_sub2 = atan2(vnyEvt_sub2[i],vnxEvt_sub2[i])/r[i].order;
	r[i].vn2d_sub2->Fill(vnxEvt_sub2[i],vnyEvt_sub2[i]);
	r[i].vn1d_sub2->Fill(sqrt(pow(vnxEvt_sub2[i],2)+pow(vnyEvt_sub2[i],2)));
	r[i].vn1dMult_sub2->Fill(sqrt(pow(vnxEvt_sub2[i],2)+pow(vnyEvt_sub2[i],2)),qncnt_sub2);
	r[i].qn_sub2+=val_sub2;
	r[i].qne_sub2+=pow(qAx_sub2,2)*qnxe_sub2/fabs(qnx_sub2) + pow(qAy_sub2,2)*qnye_sub2/fabs(qny_sub2);
	r[i].qnSub_sub2[isub]+=val_sub2;
	r[i].qnSube_sub2[isub]+=pow(qAx_sub2,2)*qnxe_sub2/fabs(qnx_sub2) + pow(qAy_sub2,2)*qnye_sub2/fabs(qny_sub2);
	r[i].wn_sub2+=qncnt_sub2;
	r[i].wne_sub2+=qncnte_sub2;

	r[i].ang_full = atan2(vnyEvt_full[i],vnxEvt_full[i])/r[i].order;
	r[i].vn2d_full->Fill(vnxEvt_full[i],vnyEvt_full[i]);
	r[i].vn1d_full->Fill(sqrt(pow(vnxEvt_full[i],2)+pow(vnyEvt_full[i],2)));
	r[i].vn1dMult_full->Fill(sqrt(pow(vnxEvt_full[i],2)+pow(vnyEvt_full[i],2)),qncnt_full);
	r[i].qn_full+=val_full;
	r[i].qne_full+=pow(qAx_full,2)*qnxe_full/fabs(qnx_full) + pow(qAy_full,2)*qnye_full/fabs(qny_full);
	r[i].qnSub_full[isub]+=val_full;
	r[i].qnSube_full[isub]+=pow(qAx_full,2)*qnxe_full/fabs(qnx_full) + pow(qAy_full,2)*qnye_full/fabs(qny_full);
	r[i].wn_full+=qncnt_full;
	r[i].wne_full+=qncnte_full;
	r[i].pt+=pt;
      }
    }

    double Ax_sub1 = qx[r[i].A1];
    double Bx_sub1 = qx[r[i].B1];
    double Cx_sub1 = qx[r[i].C1];
    double Ay_sub1 = qy[r[i].A1];
    double By_sub1 = qy[r[i].B1];
    double Cy_sub1 = qy[r[i].C1];
    double wA_sub1 = sumw[r[i].A1];
    double wB_sub1 = sumw[r[i].B1];
    double wC_sub1 = sumw[r[i].C1];
    r[i].qAB_sub1 += Ax_sub1*Bx_sub1+Ay_sub1*By_sub1;
    r[i].qAC_sub1 += Ax_sub1*Cx_sub1+Ay_sub1*Cy_sub1;
    r[i].qBC_sub1 += Bx_sub1*Cx_sub1+By_sub1*Cy_sub1;
    r[i].wAB_sub1 += wA_sub1*wB_sub1;
    r[i].wAC_sub1 += wA_sub1*wC_sub1;
    r[i].wBC_sub1 += wB_sub1*wC_sub1;

    double Ax_sub2 = qx[r[i].A2];
    double Bx_sub2 = qx[r[i].B2];
    double Cx_sub2 = qx[r[i].C2];
    double Ay_sub2 = qy[r[i].A2];
    double By_sub2 = qy[r[i].B2];
    double Cy_sub2 = qy[r[i].C2];
    double wA_sub2 = sumw[r[i].A2];
    double wB_sub2 = sumw[r[i].B2];
    double wC_sub2 = sumw[r[i].C2];
    r[i].qAB_sub2 += Ax_sub2*Bx_sub2+Ay_sub2*By_sub2;
    r[i].qAC_sub2 += Ax_sub2*Cx_sub2+Ay_sub2*Cy_sub2;
    r[i].qBC_sub2 += Bx_sub2*Cx_sub2+By_sub2*Cy_sub2;
    r[i].wAB_sub2 += wA_sub2*wB_sub2;
    r[i].wAC_sub2 += wA_sub2*wC_sub2;
    r[i].wBC_sub2 += wB_sub2*wC_sub2;

    r[i].wnA_sub1+=qncnt_sub1*wA_sub1;
    r[i].wnASub_sub1[isub]+=qncnt_sub1*wA_sub1;
    r[i].wnA_sub2+=qncnt_sub2*wA_sub2;
    r[i].wnASub_sub2[isub]+=qncnt_sub2*wA_sub2;
    r[i].wnA_full+=qncnt_sub1*wA_sub1+qncnt_sub2*wA_sub2;
    r[i].wnASub_full[isub]+=qncnt_sub1*wA_sub1+qncnt_sub2*wA_sub2;
  }
}

double Framework::GetqABC_sub1(int roi) {
  double AB = r[roi].qAB_sub1/r[roi].wAB_sub1;
  double AC = r[roi].qAC_sub1/r[roi].wAC_sub1;
  double BC = r[roi].qBC_sub1/r[roi].wBC_sub1;
  return sqrt(AB*AC/BC);
}
double Framework::GetqABC_sub2(int roi) {
  double AB = r[roi].qAB_sub2/r[roi].wAB_sub2;
  double AC = r[roi].qAC_sub2/r[roi].wAC_sub2;
  double BC = r[roi].qBC_sub2/r[roi].wBC_sub2;
  return sqrt(AB*AC/BC);
}
double Framework::GetqABC_full(int roi) {
  return (GetqABC_sub1(roi)+GetqABC_sub2(roi))/2.;
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
  double etamin1 = r[roi].minEta1;
  double etamax1 = r[roi].maxEta1;
  int ietamin1 = spec->GetYaxis()->FindBin(etamin1+0.01);
  int ietamax1 = spec->GetYaxis()->FindBin(etamax1-0.01);
  double etamin2 = r[roi].minEta2;
  double etamax2 = r[roi].maxEta2;
  int ietamin2 = spec->GetYaxis()->FindBin(etamin2+0.01);
  int ietamax2 = spec->GetYaxis()->FindBin(etamax2-0.01);

  TH1D * spec1d = (TH1D *) spec->ProjectionX(Form("spec1d_%d_%d",minCent,maxCent),ietamin1,ietamax1);
  spec1d->SetDirectory(0);
  spec1d->Add((TH1D *) spec->ProjectionX(Form("spec1d_%d_%d_sub2",minCent,maxCent),ietamin2,ietamax2));
  double bineta = 0.199;
  for(int i = 0; i<spec1d->GetNbinsX(); i++) {
    double deta = etamax1-etamin1+etamax2-etamin2;
    double dpt = spec1d->GetBinWidth(i+1);
    double binpt = spec1d->GetXaxis()->GetBinCenter(i+1);
    double eff = 1;
    if(effCorrect) {
      eff = heff->GetBinContent(heff->GetXaxis()->FindBin(bineta),heff->GetYaxis()->FindBin(binpt));
    }
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

double Framework::GetVnErrSubEvt_sub1(int roi){
  double vnav=0;
  for(int i = 0; i<10; i++) vnav+=GetVnSub_sub1(roi,i);
  vnav/=10.;
  double sig = 0;
  for(int i = 0; i<10; i++) sig+=pow(GetVnSub_sub1(roi,i)-vnav,2);
  sig = sqrt(sig/9.);
  cout<<"err: "<<vnav<<"\t"<<sig<<endl;
  return sig;
}
double Framework::GetVnErrSubEvt_sub2(int roi){
  double vnav=0;
  for(int i = 0; i<10; i++) vnav+=GetVnSub_sub2(roi,i);
  vnav/=10.;
  double sig = 0;
  for(int i = 0; i<10; i++) sig+=pow(GetVnSub_sub2(roi,i)-vnav,2);
  sig = sqrt(sig/9.);
  cout<<"err: "<<vnav<<"\t"<<sig<<endl;
  return sig;
}
double Framework::GetVnErrSubEvt_full(int roi){
  double vnav=0;
  for(int i = 0; i<10; i++) vnav+=GetVnSub_full(roi,i);
  vnav/=10.;
  double sig = 0;
  for(int i = 0; i<10; i++) sig+=pow(GetVnSub_full(roi,i)-vnav,2);
  sig = sqrt(sig/9.);
  cout<<"err: "<<vnav<<"\t"<<sig<<endl;
  return sig;
}
#endif
