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
  Framework * f = new Framework();
  for(int i = 0; i<nptbins; i++) {
    int iroi = f->SetROIRange(order, 25, 30, -0.4, 0.0, ptbins[i],ptbins[i+1]);
    if(iroi<0) return NULL;
    f->SetROIEP(iroi,HFp2,HFm2,trackmid2);
  }
  TH1D * spec = f->GetSpectra(0);
  int count = 0;
  int partcount = 0;
  for(int i = 0; i<f->GetN(); i++) {
    signal(SIGINT,intHandler);
    if(keepRunning == 0) break;
    f->AddEvent(i);
    ++count;
    ++partcount;
    if(partcount == 10000) {
      cout<<count<<endl;
      partcount = 0;
    }
  }
  cout<<"Processed "<<count<<" events"<<endl;
  double x[17];
  double xSub[17];
  double y[17];
  double ey[17];
  double eySub[17];
  for(int i = 0; i<nptbins; i++) {
    x[i] = f->GetAvPt(i);
    xSub[i] = x[i]+0.02;
    y[i] = f->GetVn(i);
    ey[i] = f->GetVnErr(i);
    eySub[i] = f->GetVnErrSubEvt(i);
    cout<<x[i]<<"\t"<<y[i]<<"\t"<<ey[i]<<"\t"<<eySub[i]<<endl;
  }
  TGraphErrors * g = new TGraphErrors(17,x,y,0,ey);
  TGraphErrors * gSub = new TGraphErrors(17,xSub,y,0,eySub);
  TCanvas * c = new TCanvas("c","c",1000,800);
  c->Divide(2);
  c->cd(1);
  TH1D * h = new TH1D("h","h",100,0,12);
  h->SetMinimum(0);
  h->SetMaximum(0.4);
  h->Draw();
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlue);
  g->SetLineColor(kBlue);
  g->Draw("p");
  gSub->SetMarkerStyle(21);
  gSub->SetMarkerColor(kRed);
  gSub->SetLineColor(kRed);
  gSub->Draw("p");
  c->cd(2);
  gPad->SetGrid();
  gPad->SetLogy();
  spec->SetMaximum(1e8);
  spec->Draw();
  TCanvas * c2 = new TCanvas("vn2d","vn2d",1200,1200);
  c2->Divide(4,4);
  for(int i = 0; i<16; i++) {
    c2->cd(i+1);
    f->Get2d(i)->Draw();
    gPad->SetGrid(1,1);
  }

  return g ;
}

#endif
