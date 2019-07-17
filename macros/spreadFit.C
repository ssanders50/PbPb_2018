#include "Math/SpecFuncMathCore.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLegend.h"

#include "/home/sanders/DAgostiniUnfold/DAgostiniUnfold.h"
#include "/home/sanders/RooUnfold/src/RooUnfoldResponse.h"
#include "/home/sanders/RooUnfold/src/RooUnfoldBayes.h"

#include "ePower.h"
#include <iostream>

TFile * tffull = new TFile("checkXY.root","read");
TFile * tfp = new TFile("checkXY_p.root","read");
TFile * tfm = new TFile("checkXY_m.root","read");
  //TFile * tf = new TFile("checkXY.root","read");

static const int order = 2;
  static const int NITER  = 6;
  static const int iter[] = {1, 2, 4, 8, 16, 32, 64, 128};
  int col[]               = {kOrange-2, kGreen+3, kCyan, kMagenta, kViolet-1, kBlue, kRed, kGray+2};

TH1D * hunfold[NITER];
TH1D * hrefold[NITER];
TH1D * hratio[NITER];
static const bool fixkn = false;
static const float fixkn_val = 0.35;
static const int nptbins = 15;
static const float ptbinsMin[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0, 2.5,  3.0,  3.5, 0.3, 0.6, 1.0};
static const float ptbinsMax[]={0.4, 0.5, 0.6,  0.8,  1.0, 1.25,  1.50,   2.0,  2.5, 3.0,  3.5,  4.0, 3.0, 3.0, 3.0};

static const int ncentbins = 13;
static const int centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,100};
static const int Order = 2;
static const int LINEWIDTH = 1;

static const double alphamin = 30;
static const double alphamax = 500;
static const double knmin = 0.1;
static const double knmax = 0.8;
static const double e0min = 0.02;
static const double e0max = 0.8;

void spreadFitRun(int minCent = 25, int maxCent=30, double minPt = 0.8, double maxPt=1.0, double vnrange=0.8, string cs="full"){
  TFile * tf = NULL;
  if(cs.find("full")!=std::string::npos) {
    tf = tffull;
  } else if (cs.find("pos")!=std::string::npos) {
    tf = tfp;
  } else if (cs.find("neg")!=std::string::npos) {
    tf = tfm;
  }
  if(tf==NULL) {
    cout<<"illegal charge state"<<endl;
    return;
  }
  ROOT::Math::MathMoreLibrary::Load();
  gaus2dFunc->SetParName(0,"vn");
  gaus2dFunc->SetParName(1,"delta");
  gaus2dFunc->SetParName(2,"scale");
  vnGausFunc->SetParName(0,"delta");
  vnGausFunc->SetParName(1,"alpha");
  vnGausFunc->SetParName(2,"kn");
  vnGausFunc->SetParName(3,"e0");
  vnGausFunc->SetParName(4,"scale");

  TH2D * hResp = (TH2D *) tf->Get(Form("%d_%d/%04.2f_%04.2f/resp2d",minCent,maxCent,minPt,maxPt));
  hResp->SetDirectory(0);

  TCanvas * c = new TCanvas(Form("c_%d_%d_%04.1f_%04.1f_%s",minCent,maxCent,minPt,maxPt, cs.data()),Form("c_%d_%d_%04.1f_%04.1f_%s",minCent,maxCent,minPt,maxPt, cs.data()),1200,1100);
  c->Divide(1,2);
  c->cd(1);
   TPad * top = (TPad * )gPad;
  top->Divide(3);
  top->cd(1);
  string sxname = Form("%d_%d/%04.2f_%04.2f/diffX",minCent,maxCent,minPt,maxPt);
  cout<<sxname<<endl;
  TH1D * diffX = (TH1D *) tf->Get(sxname.data());
  cout<<diffX<<endl;
  string syname = Form("%d_%d/%04.2f_%04.2f/diffY",minCent,maxCent,minPt,maxPt);
  TH1D * diffY = (TH1D *) tf->Get(syname.data());
  diffX->SetYTitle("Unit Normalized Counts");
  diffY->SetYTitle("Unit Normalized Counts");
  string sangcorr = Form("%d_%d/%04.2f_%04.2f/angcorr",minCent,maxCent,minPt,maxPt);
  TH2D * angcorr = (TH2D *) tf->Get(sangcorr.data());
  TF1 * fg = new TF1("fg","gaus",-0.5,0.5);
  fg->SetLineColor(kRed);
  fg->SetLineWidth(LINEWIDTH);

  cout<<"  ==========  1  ============="<<endl;
  diffX->Scale(1/diffX->Integral(1,diffX->GetNbinsX()));
  diffY->SetAxisRange(-vnrange,vnrange);
  diffY->Scale(1/diffY->Integral(1,diffY->GetNbinsX()));
  double diffm = diffX->GetMaximum();
  if(diffY->GetMaximum()>diffm) diffm=diffY->GetMaximum();
  double setm = 0.03;
  while(setm<diffm) setm+=0.01;
  diffX->SetMaximum(setm);
  diffY->SetMaximum(setm);
  diffX->SetAxisRange(-vnrange,vnrange);
  diffY->SetAxisRange(-vnrange,vnrange);
  diffX->Draw();
  gPad->SetGrid(1,1);
  diffX->Fit(fg);
  //diffX->Fit(fs,"+");
  double sigx = ((TF1 *)diffX->FindObject("fg"))->GetParameter(2);
  TLatex * lsigx = new TLatex(-0.9*vnrange,0.020*setm/0.03,Form("#sigma = %5.3f",sigx));
  TLatex * lchix = new TLatex(-0.9*vnrange,0.018*setm/0.03,Form("#chi^{2}/NDF = %5.3f",((TF1 *) diffX->FindObject("fg"))->GetChisquare()/((TF1 *) diffX->FindObject("fg"))->GetNDF()));
  lsigx->Draw();
  lchix->Draw();
  TLatex * centrange = new TLatex(-0.9*vnrange,.027*setm/0.03,Form("%d - %d%c",minCent,maxCent,'%'));
  centrange->Draw();
  TLatex * ptrange = new TLatex(-0.9*vnrange,.025*setm/0.03,Form("%04.2f < p_{T} < %04.2f GeV/c",minPt,maxPt));
  ptrange->Draw();
						     
  
  cout<<"  ==========  2  ============="<<endl;
  top->cd(2);
  diffY->Draw();
  gPad->SetGrid(1,1);
  diffY->Fit(fg);
  //diffY->Fit(fs,"+");
  double sigy = ((TF1 *)diffY->FindObject("fg"))->GetParameter(2);
  TLatex * lsigy = new TLatex(-0.9*vnrange,0.020*setm/0.03,Form("#sigma = %5.3f",sigy));
  TLatex * lchiy = new TLatex(-0.9*vnrange,0.018*setm/0.03,Form("#chi^{2}/NDF = %5.3f",((TF1 *) diffY->FindObject("fg"))->GetChisquare()/((TF1 *) diffY->FindObject("fg"))->GetNDF()));
  lsigy->Draw();
  lchiy->Draw();

  double sig = (sigx+sigy)/2.;

  
  cout<<"  ===============3 ==============="<<endl;
  top->cd(3);
  angcorr->GetXaxis()->SetRangeUser(-2,2);
  angcorr->GetYaxis()->SetRangeUser(-2,2);
  TLine * line = new TLine(-2,-2,2,2);
  line->SetLineColor(kRed);
  gPad->SetGrid(1,1);
  gPad->SetTopMargin(0.25);

  angcorr->Draw();
  line->Draw();
  
  cout<<"  ==========  4  ============="<<endl;
  c->cd(2);
  TPad * bot = (TPad *) gPad;
  bot->cd();
  
  TH1D * hSetup = new TH1D("hSetup","hSetup",1000,0,1.2);
  hSetup->SetXTitle(Form("(v_{%d}^{a} + v_{%d}^{b})/2)",order,order));
  hSetup->SetYTitle(Form("p(v_{%d})",order));
  hSetup->SetMinimum(1e-6);
  hSetup->SetMaximum(1.0);
  gPad->SetLogy();
  gPad->SetGrid(1,1);
  gPad->SetLeftMargin(0.08);
  hSetup->GetYaxis()->SetTitleOffset(0.7);
  hSetup->Draw();
  string vnname = Form("%d_%d/%04.2f_%04.2f/vn",minCent,maxCent,minPt,maxPt);
  TH1D * vn = (TH1D *) tf->Get(vnname.data());
  vn->SetYTitle("p(v_{n})");
  vn->SetLineColor(kBlue);
  vn->SetMarkerColor(kBlue);
  vn->Scale(1/vn->Integral(1,vn->GetNbinsX()));
  TH1D * hMeas = (TH1D *) vn->Clone("hMeas");
  //double scmax = vn->GetMaximum();
  //double s = 0.04;
  //while (scmax>s+0.02) s+=0.02;
  //vn->SetMaximum(s+0.04);
  vn->SetMaximum(1.0);
  vn->SetMinimum(1e-6);
  vn->Draw("same");
  cout<<"Blue circular pointer are the data points"<<endl;
  gPad->SetGrid(1,1);
  gPad->SetLogy();

  cout<<"Green dotted line is Bessel-Gaussian fit to data with width set by spreading.  2-D Gaussian spreading leads to this."<<endl;
  gaus2dFunc->SetParameter(0,.01);
  gaus2dFunc->SetLineColor(kGreen);
  gaus2dFunc->SetLineStyle(2);
  gaus2dFunc->SetLineWidth(LINEWIDTH);
  gaus2dFunc->SetRange(0,vn->GetMean()+4*sig);
  //gaus2dFunc->SetParLimits(1,0.01,0.2);
  gaus2dFunc->FixParameter(1,sig);
  gaus2dFunc->FixParameter(2,0.01);
  vn->Fit("gaus2dFunc","R");

  cout<<"Magenta line is result of elliptic power folded with Bessel-Gaussian"<<endl;
  TH1D * gaus2dfithist = (TH1D *)((TF1 *) vn->FindObject("gaus2dFunc"))->GetHistogram();
  vnGausFunc->SetParameters(sig, 58.6, .521, .210, 0.01);
  vnGausFunc->FixParameter(0,sig);
  vnGausFunc->SetParLimits(1,alphamin,alphamax);
  vnGausFunc->SetParLimits(2,knmin,knmax);
  //vnGausFunc->FixParameter(2,0.35);
  vnGausFunc->SetParLimits(3,e0min,e0max);
  vnGausFunc->FixParameter(4,0.01);
  vnGausFunc->SetLineColor(kMagenta);
  vnGausFunc->SetLineWidth(LINEWIDTH);
  vnGausFunc->SetRange(0, vn->GetMean()+3.8*sig);
  vn->Fit("vnGausFunc","+R");

  cout<<"Solid black line is elliptic power function folded with Bessel-Gaussian spread"<<endl;
  ePowerFunc->SetParameters(vnGausFunc->GetParameter(1),vnGausFunc->GetParameter(2),vnGausFunc->GetParameter(3),0.01);
  double rmax = vn->GetMean();
  while(ePowerFunc->Eval(rmax)>1e-5) rmax+=0.01;
  ePowerFunc->SetRange(0,rmax);
  ePowerFunc->Draw("same");
  
  double chi2_vnGaus = ((TF1 *) vn->FindObject("vnGausFunc"))->GetChisquare()/((TF1 *) vn->FindObject("vnGausFunc"))->GetNDF();
  double chi2_gaus = ((TF1 *) vn->FindObject("gaus2dFunc"))->GetChisquare()/((TF1 *) vn->FindObject("gaus2dFunc"))->GetNDF();


  TPaveText * pav4 = new TPaveText(0.42,0.49,0.60,0.89,"NDC");
  pav4->SetFillColor(kWhite);
  pav4->SetBorderSize(0);
  pav4->SetTextAlign(kHAlignLeft);
  pav4->AddText(Form("v_{2}         = %6.3f",((TF1 *) vn->FindObject("gaus2dFunc"))->GetParameter(0)));
  pav4->AddText(Form("#sigma          = %6.3f",((TF1 *) vn->FindObject("gaus2dFunc"))->GetParameter(1)));
  pav4->AddText("#alpha          =");
  pav4->AddText("k_{n}         =");
  pav4->AddText("#epsilon_{n}         =");
  pav4->AddText(Form("scale     = %6.3f",((TF1 *) vn->FindObject("gaus2dFunc"))->GetParameter(2)));
  pav4->AddText(Form("#chi^{2}/NDF  = %5.2f",chi2_gaus));
  pav4->SetTextFont(43);
  pav4->SetTextSize(22);
  pav4->Draw();

  TPaveText * pav5 = new TPaveText(0.62,0.5,0.77,0.90,"NDC");
  pav5->SetFillColor(kWhite);
  pav5->SetBorderSize(0);
  pav5->SetTextAlign(kHAlignLeft);

  double alpha = ((TF1 *) vn->FindObject("vnGausFunc"))->GetParameter(1);
  double ealpha = ((TF1 *) vn->FindObject("vnGausFunc"))->GetParError(1);
  double kn =   ((TF1 *) vn->FindObject("vnGausFunc"))->GetParameter(2);
  double ekn =  ((TF1 *) vn->FindObject("vnGausFunc"))->GetParError(2);
  double e0 =   ((TF1 *) vn->FindObject("vnGausFunc"))->GetParameter(3);
  double ee0 =  ((TF1 *) vn->FindObject("vnGausFunc"))->GetParError(3);

  double vnval = kn*e0;
  double vnerr = vnval*sqrt(pow(ekn/kn,2)+pow(ee0/e0,2));

  pav5->AddText(Form("%6.3f #pm %6.3f",vnval,vnerr));
  pav5->AddText(Form("%6.3f",((TF1 *) vn->FindObject("vnGausFunc"))->GetParameter(0)));

  if(fabs(alpha-alphamin)<1e-3) {
    pav5->AddText(Form("%6.3f#downarrow #pm %6.3f",alpha,ealpha));
  } else if(fabs(alpha-alphamax)<1e-3) {
    pav5->AddText(Form("%6.3f#uparrow #pm %6.3f",alpha,ealpha));
  } else {
    pav5->AddText(Form("%6.3f #pm %6.3f",alpha,ealpha));
  }
  if(fixkn) {
    pav5->AddText(Form("%6.3f",kn));
  } else {
    pav5->AddText(Form("%6.3f #pm %6.3f",kn,ekn));
  }
  pav5->AddText(Form("%6.3f #pm %6.3f",e0,ee0));
  pav5->AddText(Form("%6.3f",((TF1 *) vn->FindObject("vnGausFunc"))->GetParameter(4)));
  pav5->AddText(Form("%5.2f",chi2_vnGaus));
  pav5->SetTextFont(43);
  pav5->SetTextSize(22);
  pav5->Draw();

  RooUnfoldResponse * response = new RooUnfoldResponse(0,0,hResp,Form("response_%d_%d_%04.2f_%04.2f",minCent,maxCent,minPt,maxPt));

  DAgostiniUnfold unfold (hResp);
  //unfold.DoSystematics();

  for(int i = 0; i<NITER; i++) {
    double chi = 0;
    double ndf = 0;
    unfold.Unfold(hMeas,iter[i]);
    hunfold[i] = (TH1D *) unfold.GetHReco(Form("hreco_%d_%d_%04.2f_%04.2f_%i",minCent,maxCent,minPt,maxPt,iter[i]));
    hunfold[i]->SetLineColor(col[i]);
    hunfold[i]->SetMarkerColor(col[i]);
    //hunfold[i]->Scale(1./hunfold[i]->Integral(1,hunfold[i]->GetNbinsX()));
    hrefold[i] = unfold.refold(hunfold[i],Form("hrefold_%d_%d_%04.2f_%04.2f_%i",minCent,maxCent,minPt,maxPt,iter[i]));
    //   hrefold[i]->Sumw2();
    hratio[i] = (TH1D *) hMeas->Clone(Form("hratio_%d_%d_%04.2f_%04.2f_%i",minCent,maxCent,minPt,maxPt,iter[i]));
    hratio[i]->Divide(hrefold[i]);
    for(int j = 1; j<=hrefold[i]->GetNbinsX();j++) hrefold[i]->SetBinError(j,0);
    cout<<"iter: "<<i+1<<"  chisq/ndf = "<<hrefold[i]->Chi2Test(hMeas,"UW CHI2")<<endl;
    for(int j = 1; j<=hrefold[i]->GetNbinsX(); j++) {
      if(hMeas->GetBinContent(j)>1e-5){
	cout<<hrefold[i]->GetBinCenter(j)<<"\t"<<hrefold[i]->GetBinContent(j)<<"\t"<<hrefold[i]->GetBinError(j)<<"\t"<<hMeas->GetBinContent(j)<<"\t"<<hMeas->GetBinError(j)<<endl;
      }
    }
  }
  hunfold[NITER-1]->Draw("same");
  hrefold[NITER-1]->SetMarkerColor(kCyan);
  hrefold[NITER-1]->SetLineColor(kCyan);
  hrefold[NITER-1]->Draw("same");
  cout<<"Fit Unfolded Distribution"<<endl;
  TF1 * newePowerFunc = new TF1("newePowerFunc",ePowerWrap,0,1,4);
  newePowerFunc->SetParameters(alpha,kn,e0,0.01);
  newePowerFunc->FixParameter(3,0.01);
  newePowerFunc->SetParLimits(0,alphamin,alphamax);
  newePowerFunc->SetParLimits(1,knmin,knmax);
  newePowerFunc->SetParLimits(2,e0min,e0max);
  newePowerFunc->SetLineColor(kCyan+2);
  newePowerFunc->SetLineWidth(LINEWIDTH);
  newePowerFunc->SetRange(0,rmax);
  hunfold[NITER-1]->Fit(newePowerFunc,"R");

  TPaveText * pav6 = new TPaveText(0.79,0.5,0.94,0.90,"NDC");
  pav6->SetFillColor(kWhite);
  pav6->SetBorderSize(0);
  pav6->SetTextAlign(kHAlignLeft);
  double chi2_ep = newePowerFunc->GetChisquare()/newePowerFunc->GetNDF();

  double nalpha = newePowerFunc->GetParameter(0);
  double nealpha = newePowerFunc->GetParError(0);
  double nkn =     newePowerFunc->GetParameter(1);
  double nekn =    newePowerFunc->GetParError(1);
  double ne0 =     newePowerFunc->GetParameter(2);
  double nee0 =    newePowerFunc->GetParError(2);

  double nvnval = kn*e0;
  double nvnerr = nvnval*sqrt(pow(nekn/nkn,2)+pow(nee0/ne0,2));

  pav6->AddText(Form("%6.3f #pm %6.3f",nvnval,nvnerr));
  pav6->AddText("  ");

  if(fabs(nalpha-alphamin)<1e-3) {
    pav6->AddText(Form("%6.3f#downarrow #pm %6.3f",nalpha,nealpha));
  } else if(fabs(nalpha-alphamax)<1e-3) {
    pav6->AddText(Form("%6.3f#uparrow #pm %6.3f",nalpha,nealpha));
  } else {
    pav6->AddText(Form("%6.3f #pm %6.3f",nalpha,nealpha));
  }
  if(fixkn) {
    pav6->AddText(Form("%6.3f",nkn));
  } else {
    pav6->AddText(Form("%6.3f #pm %6.3f",nkn,nekn));
  }
  pav6->AddText(Form("%6.3f #pm %6.3f",e0,ee0));
  pav6->AddText(Form("%6.3f",newePowerFunc->GetParameter(3)));
  pav6->AddText(Form("%5.2f",chi2_ep));
  pav6->SetTextFont(43);
  pav6->SetTextSize(22);
  pav6->SetTextColor(kCyan+2);
  pav6->Draw();


  TLegend * leg = new TLegend(0.68,0.2,0.9,0.45);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
  leg->SetTextSize(14);
  leg->AddEntry(vn,"v_{n}^{obs}","lp");
  leg->AddEntry(vnGausFunc,"Elliptic Power w/ spread. #sigma_{spread} fixed.","l");
  leg->AddEntry(gaus2dfithist,"2D Gaussian.  #sigma = #sigma_{spread}.","l");
  leg->AddEntry(ePowerFunc,"Elliptic Power w/o spread.","l");
  leg->AddEntry(hunfold[NITER-1],"Unfolded distribution","lp");
  leg->AddEntry(newePowerFunc,"Elliptic Power fit to Unfolded","l");
  leg->Draw();
  string canname = "figures/spreadFit/"+to_string(minCent)+"_"+to_string(maxCent)+"/";
  c->Print(Form("%s/pt_%4.2f_%4.2f_%s.pdf",canname.data(),minPt,maxPt,cs.data()),"pdf");
	     
}

void spreadFit(){
  //for(int k = 0; k<ncentbins; k++) {
  //  for(int i = 0; i<nptbins; i++) {
  //   spreadFitRun(centbins[k],centbins[k+1],ptbinsMin[i],ptbinsMax[i],0.8,"full");
  //   spreadFitRun(centbins[k],centbins[k+1],ptbinsMin[i],ptbinsMax[i],0.8,"pos");
  //   spreadFitRun(centbins[k],centbins[k+1],ptbinsMin[i],ptbinsMax[i],0.8,"neg");
  // }
  //}
   spreadFitRun(20,25,0.3,3.0,0.8,"full");
   //spreadFitRun(20,25,0.6,3.0,0.8,"pos");
   //spreadFitRun(20,25,0.6,3.0,0.8,"neg");

}
