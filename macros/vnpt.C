#ifndef VNPT
#define VNPT
#include <signal.h>
#include "Framework.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "../CMSSW_10_3_2/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;
using namespace std;
static volatile int keepRunning = 1;
void intHandler(int dummy){
  keepRunning = 0;
}
static const int nptbins = 17;
static const float ptbins[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0,
			      2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  7.0, 8.0, 10.};

TGraphErrors * vnpt(){
  int order = 2;
  Framework * fm = new Framework();
  Framework * fp = new Framework();
  for(int i = 0; i<nptbins; i++) {
    int iroim = fm->SetROIRange(order, 25, 30, -0.4, 0.0, ptbins[i],ptbins[i+1]);
    if(iroim<0) return NULL;
    fm->SetROIEP(iroim,HFp2,HFm2,trackmid2);
    int iroip = fp->SetROIRange(order, 25, 30, 0.0, 0.4, ptbins[i],ptbins[i+1]);
    if(iroip<0) return NULL;
    fp->SetROIEP(iroip,HFm2,HFp2,trackmid2);
  }
  TH1D * specm = fm->GetSpectra(0);
  TH1D * specp = fp->GetSpectra(0);
  int count = 0;
  int partcount = 0;
  for(int i = 0; i<fm->GetN(); i++) {
    signal(SIGINT,intHandler);
    if(keepRunning == 0) break;
    fm->AddEvent(i);
    double fmx = fm->GetVnxEvt();
    double fmy = fm->GetVnyEvt();
    fp->AddEvent(i);
    double fpx = fp->GetVnxEvt();
    double fpy = fp->GetVnyEvt();

    ++count;
    ++partcount;
    if(partcount == 10000) {
      cout<<count<<endl;
      partcount = 0;
    }
  }
  cout<<"Processed "<<count<<" events"<<endl;
  //
  double xm[17];
  double xmSub[17];
  double ym[17];
  double eym[17];
  double eymSub[17];
  for(int i = 0; i<nptbins; i++) {
    xm[i] = fm->GetAvPt(i);
    xmSub[i] = xm[i]+0.02;
    ym[i] = fm->GetVn(i);
    eym[i] = fm->GetVnErr(i);
    eymSub[i] = fm->GetVnErrSubEvt(i);
    cout<<xm[i]<<"\t"<<ym[i]<<"\t"<<eym[i]<<"\t"<<eymSub[i]<<endl;
  }
  TGraphErrors * gm = new TGraphErrors(17,xm,ym,0,eym);
  TGraphErrors * gmSub = new TGraphErrors(17,xmSub,ym,0,eymSub);
  TCanvas * c = new TCanvas("c","c",1000,800);
  c->Divide(2);
  c->cd(1);
  TH1D * h = new TH1D("h","h",100,0,12);
  h->SetMinimum(0);
  h->SetMaximum(0.4);
  h->Draw();
  gm->SetMarkerStyle(20);
  gm->SetMarkerColor(kBlue);
  gm->SetLineColor(kBlue);
  gm->Draw("p");
  gmSub->SetMarkerStyle(21);
  gmSub->SetMarkerColor(kRed);
  gmSub->SetLineColor(kRed);
  gmSub->Draw("p");
  //
  double xp[17];
  double xpSub[17];
  double yp[17];
  double eyp[17];
  double eypSub[17];
  for(int i = 0; i<nptbins; i++) {
    xp[i] = fp->GetAvPt(i);
    xpSub[i] = xp[i]+0.02;
    yp[i] = fp->GetVn(i);
    eyp[i] = fp->GetVnErr(i);
    eypSub[i] = fp->GetVnErrSubEvt(i);
    cout<<xp[i]<<"\t"<<yp[i]<<"\t"<<eyp[i]<<"\t"<<eypSub[i]<<endl;
  }
  TGraphErrors * gp = new TGraphErrors(17,xp,yp,0,eyp);
  TGraphErrors * gpSub = new TGraphErrors(17,xpSub,yp,0,eypSub);
  //TCanvas * c = new TCanvas("c","c",1000,800);
  //c->Divide(2);
  //c->cd(1);
  //TH1D * h = new TH1D("h","h",100,0,12);
  //h->SetMinimum(0);
  //h->SetMaximum(0.4);
  //h->Draw();
  gp->SetMarkerStyle(20);
  gp->SetMarkerColor(kGreen);
  gp->SetLineColor(kGreen);
  gp->Draw("p");
  gpSub->SetMarkerStyle(21);
  gpSub->SetMarkerColor(kCyan);
  gpSub->SetLineColor(kCyan);
  gpSub->Draw("p");

  //
  c->cd(2);
  gPad->SetGrid();
  gPad->SetLogy();
  specm->SetMaximum(1e8);
  specm->Draw();
  specp->SetLineColor(kRed);
  specp->SetMarkerColor(kRed);
  specp->Draw("same");
  TCanvas * c2 = new TCanvas("vn2d","vn2d",1200,1200);
  c2->Divide(4,4);
  for(int i = 0; i<16; i++) {
    c2->cd(i+1);
    fm->Get2d(i)->Draw();
    gPad->SetGrid(1,1);
  }

  return gm ;
}

#endif
