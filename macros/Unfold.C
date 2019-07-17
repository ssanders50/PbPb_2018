#include <iostream>
using namespace std;
#include "TSystem.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2.h"
#include "TFile.h"
#include "cmath"

#include "/home/sanders/DAgostiniUnfold/DAgostiniUnfold.h"
#include "/home/sanders/RooUnfold/src/RooUnfoldResponse.h"
#include "/home/sanders/RooUnfold/src/RooUnfoldBayes.h"

static const int order = 2;
static const int NITER  = 6;
static const int iter[] = {1, 2, 4, 8, 16, 32, 64, 128};
int col[]               = {kOrange-2, kGreen+3, kCyan, kMagenta, kViolet-1, kBlue, kRed, kGray+2};

TH1D * hunfold[NITER];
TH1D * hrefold[NITER];
TH1D * hratio[NITER];

TFile * fin = new TFile("checkXY.root");

const Double_t cutdummy= -99999.0;

// Double_t smear (Double_t xt)
// {
//   Double_t x= gRandom->Rndm();
//   Double_t xsmear= gRandom->Gaus(0,2.);     // bias and smear
//   return xt+xsmear;
// }

void Unfold(int centMin=15, int centMax=20, double ptMin=0.3, double ptMax=3.0){
  TH2D * hResp = (TH2D *) fin->Get(Form("%d_%d/%04.2f_%04.2f/resp2d",centMin,centMax,ptMin,ptMax));
  hResp->SetDirectory(0);
  TH1D * hMeas = (TH1D *) fin->Get(Form("%d_%d/%04.2f_%04.2f/vn",centMin,centMax,ptMin,ptMax));
  hMeas->Sumw2();
  hMeas->SetDirectory(0);

  DAgostiniUnfold unfolder(hResp);
  //unfolder.DebugMessages();
  //unfolder.SetSmoothIter();
  //unfolder.DoSystematics();
  unfolder.Unfold(hMeas,16);
  TH1D * hReco = (TH1D*) unfolder.GetHReco("reco");
  TH1D * hRefold = (TH1D *) unfolder.refold(hReco,"refold");
  hReco->SetLineColor(kRed);
  hReco->SetMarkerColor(kRed);
  hReco->Draw();
  hMeas->SetLineColor(kBlue);
  hMeas->SetMarkerColor(kBlue);
  hMeas->Draw("same");
  hRefold->SetLineColor(kGreen);
  hRefold->SetMarkerColor(kGreen);
  hRefold->Draw("same");
  for(int i = 1; i<=hRefold->GetNbinsX(); i++) if(hMeas->GetBinContent(i)==0) hRefold->SetBinContent(i,0);
  double chi2 = hRefold->Chi2Test(hMeas,"CHI2/NDF");
  cout<<"chi2: "<<chi2<<endl;
}
