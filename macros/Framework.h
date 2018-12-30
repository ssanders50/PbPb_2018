#ifndef FRAMEWORK
#define FRAMEWORK
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TRandom.h"
#include "../CMSSW_10_3_2/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;
using namespace std;
class Framework{
public:
  Framework(string fileName="data.root");
  int GetN(){return maxevents;}
  void GetEntry(int i){tr->GetEntry(i);}
  int SetROIRange(int order, int minCent, int maxCent, double minEta, double maxEta, double minPt, double maxPt);
  void ShowAllROIRanges();
  void SetROIEP(int roi, int EPA, int EPB, int EPC, int EPA2=-1, int EPB2=-1, int EPC2 = -1);
  void AddEvent(int evt);
  double GetAvPt(int roi){return r[roi].pt/r[roi].wn;}
  double GetqnA(int roi){return r[roi].qn/r[roi].wnA;}
  double GetqABC(int roi);
  double GetSpectraRaw(int roi);
private:
  TRandom * ran;
  int maxevents;
  TFile * tf;
  TTree * tr;
  double centval;
  int ntrkval;
  int Noff;
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
  struct range {
    int order;
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
    double pt=0;
    double wn=0;
    double wnA = 0;
    double qAB = 0;
    double wAB = 0;
    double qAC = 0;
    double wAC = 0;
    double qBC = 0;
    double wBC = 0;
    double qnSub[10]={0,0,0,0,0,0,0,0,0,0};
    double wnSub[10]={0,0,0,0,0,0,0,0,0,0};
  } r[20];
  int nrange = 0;
};

Framework::Framework(string fileName){
  ran = new TRandom();
  tf = new TFile(fileName.data(),"read");
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
  tr->SetBranchAddress("qxtrk1",  &qxtrk[0]);
  tr->SetBranchAddress("qytrk1",  &qytrk[0]);
  tr->SetBranchAddress("qxtrk2",  &qxtrk[1]);
  tr->SetBranchAddress("qytrk2",  &qytrk[1]);
  tr->SetBranchAddress("qxtrk3",  &qxtrk[2]);
  tr->SetBranchAddress("qytrk3",  &qytrk[2]);
  tr->SetBranchAddress("qxtrk4",  &qxtrk[3]);
  tr->SetBranchAddress("qytrk4",  &qytrk[3]);
  tr->SetBranchAddress("qxtrk5",  &qxtrk[4]);
  tr->SetBranchAddress("qytrk5",  &qytrk[4]);
  tr->SetBranchAddress("qxtrk6",  &qxtrk[5]);
  tr->SetBranchAddress("qytrk6",  &qytrk[5]);
  tr->SetBranchAddress("qxtrk7",  &qxtrk[6]);
  tr->SetBranchAddress("qytrk7",  &qytrk[6]);
  tr->SetBranchAddress("qcnt",    &qcnt);
  tr->SetBranchAddress("avpt",    &avpt);
  
}
int Framework::SetROIRange(int order, int minCent, int maxCent, double minEta, double maxEta, double minPt, double maxPt) {
  if(maxevents==0) return -1;
  GetEntry(0);
  r[nrange].order = order;
  r[nrange].minCent = minCent;
  r[nrange].maxCent = maxCent;
  r[nrange].minEtaBin = qxtrk[0]->GetYaxis()->FindBin(minEta);
  r[nrange].maxEtaBin = qxtrk[0]->GetYaxis()->FindBin(maxEta-0.01);
  r[nrange].minPtBin = qxtrk[0]->GetXaxis()->FindBin(minPt);
  r[nrange].maxPtBin = qxtrk[0]->GetXaxis()->FindBin(maxPt-0.01);
  r[nrange].minEta = qxtrk[0]->GetYaxis()->GetBinLowEdge(r[nrange].minEtaBin);
  r[nrange].maxEta = qxtrk[0]->GetYaxis()->GetBinLowEdge(r[nrange].maxEtaBin)+qxtrk[0]->GetYaxis()->GetBinWidth(r[nrange].maxEtaBin);
  r[nrange].minPt = qxtrk[0]->GetXaxis()->GetBinLowEdge(r[nrange].minPtBin);
  r[nrange].maxPt = qxtrk[0]->GetXaxis()->GetBinLowEdge(r[nrange].maxPtBin)+qxtrk[0]->GetXaxis()->GetBinWidth(r[nrange].maxPtBin);
  ++nrange;
  return nrange-1;
} 

void Framework::ShowAllROIRanges(){
  std::cout<<"indx\tminCent\tmaxCent\tminEta\tmaxEta\tminPt\tmaxPt\tminEtaBin\tmaxEtaBin\tminPtBin\tmaxPtBin"<<std::endl;
  for(int i = 0; i<nrange; i++) {
    std::cout<<i<<"\t"<<r[i].minCent<<"\t"<<r[i].maxCent<<"\t"<<r[i].minEta<<"\t"<<r[i].maxEta<<"\t"<<r[i].minPt<<"\t"<<r[i].maxPt<<"\t"<<r[i].minEtaBin<<"\t"<<r[i].maxEtaBin<<"\t"<<r[i].minPtBin<<"\t"<<r[i].maxPtBin<<std::endl;
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
  for(int i = 0; i<nrange; i++) {
    if(r[i].A<0) {
      cout<<"Event planes have not yet been initialized for ROI # "<<i<<endl;
    }
    if(centval<(double) r[i].minCent || centval>=(double)r[i].maxCent) continue;
    double qAx = qx[r[i].A];
    double qAy = qy[r[i].A];
    double qnx = 0;
    double qny = 0;
    double qncnt = 0;
    double pt = 0;
    for(int j = r[i].minEtaBin; j<= r[i].maxEtaBin; j++) {
      for(int k = r[i].minPtBin; k<=r[i].maxPtBin; k++) {
	qnx+=qxtrk[r[i].order-1]->GetBinContent(k,j);
	qny+=qytrk[r[i].order-1]->GetBinContent(k,j);
	qncnt+=qcnt->GetBinContent(k,j);
	pt+=avpt->GetBinContent(k,j);
      }
    }
    int isub = ran->Uniform(0,9.9999);
    double val = qAx*qnx+qAy*qny;
    r[i].qn+=val;
    r[i].qnSub[isub]+=val;
    r[i].wn+=qncnt;
    r[i].pt+=pt;
    r[i].wnSub[isub]+=qncnt;
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
  }

}
double Framework::GetqABC(int roi) {
  double AB = r[roi].qAB/r[roi].wAB;
  double AC = r[roi].qAC/r[roi].wAC;
  double BC = r[roi].qBC/r[roi].wBC;
  return sqrt(AB*AC/BC);
}

#endif
