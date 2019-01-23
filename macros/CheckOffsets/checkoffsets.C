#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TFile.h"
#include <iostream>
#include "HiEvtPlaneList.h"
FILE * dlist;
using namespace hi;
using namespace std;
TH1D * runcnt;
TH1D * qxOrig[NumEPNames];
TH1D * qyOrig[NumEPNames];
TH1D * qOrig[NumEPNames];
TH1D * qxFlat[NumEPNames];
TH1D * qyFlat[NumEPNames];
TH1D * qFlat[NumEPNames];
void checkoffsets(){
  float Cent;
  float Vtx;
  unsigned int Run;
  float EPAngs[NumEPNames];
  float EPOrig[NumEPNames];
  double qx[NumEPNames];
  double qy[NumEPNames];
  int minRun = 326400;
  int maxRun = 327500;
  int nbins = maxRun-minRun;
  cout<<"nbins: "<<nbins<<endl;
  runcnt = new TH1D("runcnt","runcnt",nbins,minRun,maxRun);
  for(int i = 0; i< NumEPNames; i++) {
    string xnameOrig = "qxOrig_"+EPNames[i];
    string ynameOrig = "qyOrig_"+EPNames[i];
    string nameOrig  = "qOrig_"+EPNames[i];
    string xnameFlat = "qxFlat_"+EPNames[i];
    string ynameFlat = "qyFlat_"+EPNames[i];
    string nameFlat  = "qFlat_"+EPNames[i];
    qxOrig[i]=new TH1D(xnameOrig.data(),xnameOrig.data(),nbins,minRun,maxRun);
    qyOrig[i]=new TH1D(ynameOrig.data(),ynameOrig.data(),nbins,minRun,maxRun);
    qOrig[i]=new TH1D(nameOrig.data(),nameOrig.data(),nbins,minRun,maxRun);
    qxFlat[i]=new TH1D(xnameFlat.data(),xnameFlat.data(),nbins,minRun,maxRun);
    qyFlat[i]=new TH1D(ynameFlat.data(),ynameFlat.data(),nbins,minRun,maxRun);
    qFlat[i]=new TH1D(nameFlat.data(),nameFlat.data(),nbins,minRun,maxRun);
  }
  dlist=fopen("datafiles.lis","r");
  char buf[80];

  while(fgets(buf,80,dlist)!=NULL) {
    buf[strlen(buf)-1]=0;
    TFile * f = new TFile(buf,"read");
    TTree * tree = (TTree *) f->Get("EPtree");
    tree->SetBranchAddress("Cent",    &Cent);
    tree->SetBranchAddress("Vtx",     &Vtx);
    tree->SetBranchAddress("Run",     &Run);
    tree->SetBranchAddress("EPAngs", &EPAngs);
    tree->SetBranchAddress("EPOrig", &EPOrig);
    cout<<tree->GetEntries()<<endl;
    for(int i = 0; i<tree->GetEntries(); i++){
      tree->GetEntry(i);
      if(Cent>80) continue;
      if(fabs(Vtx)>15) continue;
      runcnt->Fill(Run);
      for(int j = 0; j< NumEPNames; j++) {
	qxOrig[j]->Fill(Run,TMath::Cos(EPOrder[j]*EPOrig[j]));
	qyOrig[j]->Fill(Run,TMath::Sin(EPOrder[j]*EPOrig[j]));
	qxFlat[j]->Fill(Run,TMath::Cos(EPOrder[j]*EPAngs[j]));
	qyFlat[j]->Fill(Run,TMath::Sin(EPOrder[j]*EPAngs[j]));
      }
    }
    f->Close();
  }
  TFile * resfile = new TFile("results.root","recreate");
  resfile->cd();
  runcnt->Write();
  for(int i = 0; i<NumEPNames; i++) {
    qxOrig[i]->Divide(runcnt);
    qyOrig[i]->Divide(runcnt);
    qxFlat[i]->Divide(runcnt);
    qyFlat[i]->Divide(runcnt);
    for(int j = 1; j<= qxOrig[i]->GetNbinsX(); j++){
      double v = pow(qxOrig[i]->GetBinContent(j),2)+pow(qyOrig[i]->GetBinContent(j),2);
      if(v>0) qOrig[i]->SetBinContent(j,TMath::Sqrt(v));
      v = pow(qxFlat[i]->GetBinContent(j),2)+pow(qyFlat[i]->GetBinContent(j),2);
      if(v>0) qFlat[i]->SetBinContent(j,TMath::Sqrt(v));
    }
    qxOrig[i]->Write();
    qyOrig[i]->Write();
    qOrig[i]->Write();
    qxFlat[i]->Write();
    qyFlat[i]->Write();
    qFlat[i]->Write();
  }
  resfile->Close();
}
