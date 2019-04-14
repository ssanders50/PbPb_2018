#ifndef VNPT
#define VNPT
#include <signal.h>
#include "Framework.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveText.h"
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

TH1D * diffy[17];
TH1D * diffx[17];

TGraphErrors * vnpt(int mincent=5,int maxcent = 10, double etamin = 0, double etamax = 0.8){
  int order = 2;
  for(int i = 0; i<nptbins; i++) {
    diffx[i]=new TH1D(Form("diffx%d",i),Form("diffx%d",i),500,-1,1);
    diffy[i]=new TH1D(Form("diffy%d",i),Form("diffy%d",i),500,-1,1);
    diffx[i]->SetXTitle("(v_{2,x}^{obs,hfp}-v_{2,x}^{obs,hfm})/2");
    diffx[i]->SetYTitle("Events");
    diffy[i]->SetXTitle("(v_{2,y}^{obs,hfp}-v_{2,y}^{obs,hfm})/2");
    diffy[i]->SetYTitle("Events");

  }
  Framework * fm = new Framework();
  Framework * fp = new Framework();
  int iroim = 0;
  int iroip = 0;
  for(int i = 0; i<nptbins; i++) {
    iroim = fm->SetROIRange(order, mincent, maxcent, -etamax, -etamin, ptbins[i],ptbins[i+1]);
    if(iroim<0) return NULL;
    fm->SetROIEP(iroim,HFp2,HFm2,trackmid2);
    iroip = fp->SetROIRange(order, mincent, maxcent, etamin, etamax, ptbins[i],ptbins[i+1]);
    if(iroip<0) return NULL;
    fp->SetROIEP(iroip,HFm2,HFp2,trackmid2);
  }
  TH1D * specm = fm->GetSpectra(0);
  TH1D * specp = fp->GetSpectra(0);
  int count = 0;
  int partcount = 0;
 rep:  for(int i = 0; i<fm->GetN(); i++) {
    signal(SIGINT,intHandler);
    if(keepRunning == 0) break;
    fm->AddEvent(i);
    fp->AddEvent(i);
    for (int j = 0; j<iroim; j++){
      double fmx = fm->GetVnxEvt(j);
      double fmy = fm->GetVnyEvt(j);
      double fpx = fp->GetVnxEvt(j);
      double fpy = fp->GetVnyEvt(j);
      if(fmx>-2&&fmy>-2&&fpx>-2&&fpy>-2) {
	diffx[j]->Fill((fpx-fmx)/2.);
	diffy[j]->Fill((fpy-fmy)/2.);
      }
    }
    ++count;
    ++partcount;
    if(partcount == 50000) {
      cout<<count<<endl;
      partcount = 0;
    }
  }
  if(fm->AddFile()&&keepRunning == 1) goto rep;
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
  TCanvas * c2 = new TCanvas("vnm2d","vnm2d",1200,1200);
  c2->Divide(4,4);
  for(int i = 0; i<16; i++) {
    c2->cd(i+1);
    fm->Get2d(i)->Draw();
    gPad->SetGrid(1,1);
    TLatex * l = new TLatex(-1.2,-1.2,Form("%3.1f<p_{T}<%3.1f GeV/c",ptbins[i],ptbins[i+1]));
    l->Draw();
    if(i==0) {
      TLatex * lt = new TLatex(-1.2,1.2,Form("%3.1f < #eta < %3.1f; %d - %d%c",-etamax, etamin, mincent,maxcent,'%'));
      lt->Draw();
    }
  }

  TCanvas * c3 = new TCanvas("vnp2d","vnp2d",1200,1200);
  c3->Divide(4,4);
  for(int i = 0; i<16; i++) {
    c3->cd(i+1);
    fp->Get2d(i)->Draw();
    gPad->SetGrid(1,1);
    TLatex * l = new TLatex(-1.2,-1.2,Form("%3.1f<p_{T}<%3.1f GeV/c",ptbins[i],ptbins[i+1]));
    l->Draw();
    if(i==0) {
      TLatex * lt = new TLatex(-1.2,1.2,Form("%3.1f < #eta < %3.1f; %d - %d%c",etamin, etamax, mincent,maxcent,'%'));
      lt->Draw();
    }
  }

  TCanvas * c4 = new TCanvas("diffx","diffx",1200,1200);
  c4->Divide(4,4);
  for(int i = 0; i<16; i++) {
    c4->cd(i+1);
    diffx[i]->Draw();
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    TPaveText * pt = new TPaveText(0.2,0.2,0.6,0.27,"NDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(kWhite);
    pt->AddText(Form("%3.1f<p_{T}<%3.1f GeV/c",ptbins[i],ptbins[i+1]));
    pt->Draw();
    TPaveText * pt2 = new TPaveText(0.2,0.85,0.6,0.93,"NDC");
    pt2->SetBorderSize(0);
    pt2->SetFillColor(kWhite);
    pt2->AddText(Form("%3.1f < #eta < %3.1f; %d - %d%c",etamin, etamax, mincent,maxcent,'%'));
    pt2->Draw();
  }

  TCanvas * c5 = new TCanvas("diffy","diffy",1200,1200);
  c5->Divide(4,4);
  for(int i = 0; i<16; i++) {
    c5->cd(i+1);
    diffy[i]->Draw();
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    TPaveText * pt = new TPaveText(0.2,0.2,0.6,0.27,"NDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(kWhite);
    pt->AddText(Form("%3.1f<p_{T}<%3.1f GeV/c",ptbins[i],ptbins[i+1]));
    pt->Draw();
    TPaveText * pt2 = new TPaveText(0.2,0.85,0.6,0.93,"NDC");
    pt2->SetBorderSize(0);
    pt2->SetFillColor(kWhite);
    pt2->AddText(Form("%3.1f < #eta < %3.1f; %d - %d%c",etamin, etamax, mincent,maxcent,'%'));
    pt2->Draw();
  }
  return gm ;
}

#endif
