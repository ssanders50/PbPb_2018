TFile * f = new TFile("results.root","read");
int iovs[51]={326499, 326500, 326501, 326502, 326504, 326518, 326521, 326545, 
	    326547, 326550, 326574, 326575, 326576, 326580, 326585, 326586,
	    326587, 326617, 326620, 326627, 326635, 326790, 326822, 326854, 
	    326855, 326856, 326859, 326883, 326897, 326941, 326961, 326965, 
	    326998, 327004, 327078, 327123, 327125, 327147, 327148, 327173, 
	    327174, 327211, 327219, 327220, 327230, 327264, 327327, 327400, 
	    327403, 327424, 328000};

void plot_trackmid2(string epname = "trackmid2"){
  string sqxOrig = "qxOrig_"+epname;
  string sqyOrig = "qyOrig_"+epname;
  string sqxFlat = "qxFlat_"+epname;
  string sqyFlat = "qyFlat_"+epname;
  TH1D * qxOrig = (TH1D *) f->Get(sqxOrig.data());
  TH1D * qyOrig = (TH1D *) f->Get(sqyOrig.data());
  TH1D * rotOrig = (TH1D *) qyOrig->Clone("rotOrig");
  double avrot = 0;
  int nrot = 0;
  for(int i = 1; i<=qyOrig->GetNbinsX(); i++) {
    double y = qyOrig->GetBinContent(i);
    double ey = qyOrig->GetBinError(i);
    double x = qxOrig->GetBinContent(i);
    double ex = qxOrig->GetBinError(i);
    if(x*x+y*y > 0) {
      double ang = TMath::ATan2(y,x)/2.;
      double terr[4];
      terr[0] = fabs(TMath::ATan2(y+ey,x+ex)/2. - ang);
      terr[1] = fabs(TMath::ATan2(y-ey,x+ex)/2. - ang);
      terr[2] = fabs(TMath::ATan2(y+ey,x-ex)/2. - ang);
      terr[3] = fabs(TMath::ATan2(y-ey,x-ex)/2. - ang);
      double err = 0;
      for(int j = 0; j<4; j++) if(terr[j]>err) err=terr[j];

      rotOrig->SetBinContent(i,ang);
      rotOrig->SetBinError(i,err);
      avrot+=ang;
      ++nrot;
    }
  }
  avrot/=nrot;
  for(int i = 1; i<qyOrig->GetNbinsX(); i++) { 
    rotOrig->SetBinContent(i,(rotOrig->GetBinContent(i)-avrot)*180./3.1415);
    rotOrig->SetBinError(i, rotOrig->GetBinError(i)*180./3.1415);
  }
  TH1D * qxFlat = (TH1D *) f->Get(sqxFlat.data());
  TH1D * qyFlat = (TH1D *) f->Get(sqyFlat.data());
  TH1D * rotFlat = (TH1D *) qyFlat->Clone("rotFlat");
  avrot = 0;
  nrot = 0;
  for(int i = 1; i<=qyFlat->GetNbinsX(); i++) {
    double y = qyFlat->GetBinContent(i);
    double ey = qyFlat->GetBinError(i);
    double x = qxFlat->GetBinContent(i);
    double ex = qxFlat->GetBinError(i);
    if(x*x+y*y > 0) {
      double ang = TMath::ATan2(y+.1,x+0.1)/2.;
      double terr[4];
      terr[0] = fabs(TMath::ATan2(y+ey+.1,x+ex+0.1)/2. - ang);
      terr[1] = fabs(TMath::ATan2(y-ey+.1,x+ex+0.1)/2. - ang);
      terr[2] = fabs(TMath::ATan2(y+ey+.1,x-ex+0.1)/2. - ang);
      terr[3] = fabs(TMath::ATan2(y-ey+.1,x-ex+0.1)/2. - ang);
      double err = 0;
      for(int j = 0; j<4; j++) if(terr[j]>err) err=terr[j];

      rotFlat->SetBinContent(i,ang);
      rotFlat->SetBinError(i,err);
      avrot+=ang;
      ++nrot;
    }
  }
  avrot/=nrot;
  for(int i = 1; i<qyFlat->GetNbinsX(); i++) { 
    rotFlat->SetBinContent(i,(rotFlat->GetBinContent(i)-avrot)*180./3.1415);
    rotFlat->SetBinError(i, rotFlat->GetBinError(i)*180./3.1415);
  }
  double msize = 1.1;
  qxOrig->SetMarkerColor(kRed);
  qxOrig->SetLineColor(kRed);
  qxOrig->SetMarkerStyle(24);
  qxOrig->SetMarkerSize(msize);

  qyOrig->SetMarkerColor(kBlue);
  qyOrig->SetLineColor(kBlue);
  qyOrig->SetMarkerStyle(25);
  qyOrig->SetMarkerSize(msize);

  qxFlat->SetMarkerColor(kRed);
  qxFlat->SetLineColor(kRed);
  qxFlat->SetMarkerStyle(24);
  qxFlat->SetMarkerSize(msize);

  qyFlat->SetMarkerColor(kBlue);
  qyFlat->SetLineColor(kBlue);
  qyFlat->SetMarkerStyle(25);
  qyFlat->SetMarkerSize(msize);

  string cname = "c_"+epname;
  TCanvas * c = new TCanvas(cname.data(),cname.data(),1000,800);
  c->Divide(1,2);
  c->cd(1);
  gPad->SetLeftMargin(0.1);
  double range = 0.6;
  qxOrig->SetMinimum(-range);
  qxOrig->SetMaximum(range);
  qxOrig->SetXTitle("Run #");
  qxOrig->SetYTitle("Offset (before recentering)");
  qxOrig->GetYaxis()->SetTitleOffset(0.8);
  qxOrig->Draw();
  qyOrig->Draw("same");
  TLatex * lep = new TLatex(326500,0.15*range/0.20,epname.data());
  lep->SetTextFont(43);
  lep->SetTextSize(32);
  lep->Draw();
  TLine * l[51];
  for(int i = 0; i<51; i++) {
    l[i]=new TLine(iovs[i],range,iovs[i],+0.25*range/0.20);
    l[i]->SetLineStyle(1);
    l[i]->Draw();
  }
  TLatex * ll = new TLatex(326400,0.205*range/0.20,"iovs:");
  ll->SetTextFont(43);
  ll->SetTextSize(18);
  ll->Draw();
  c->cd(2);
  gPad->SetLeftMargin(0.1);
  qxFlat->SetMinimum(-0.2);
  qxFlat->SetMaximum(0.2);
  qxFlat->SetXTitle("Run #");
  qxFlat->SetYTitle("Offset (after recentering)");
  qxFlat->GetYaxis()->SetTitleOffset(0.8);
  qxFlat->Draw();
  qyFlat->Draw("same");
  TLegend * leg = new TLegend(0.2,0.2,0.4,0.4);
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
  leg->SetTextSize(24);
  leg->AddEntry(qxFlat,"<cos(2 #Psi_{2})>","lp");
  leg->AddEntry(qyFlat,"<sin(2 #Psi_{2})>","lp");
  leg->Draw();
  for(int i = 0; i<51; i++) {
    l[i]=new TLine(iovs[i],0.2,iovs[i],+0.25);
    l[i]->SetLineStyle(1);
    l[i]->Draw();
  }
  TLatex * ll2 = new TLatex(326400,0.205,"iovs:");
  ll2->SetTextFont(43);
  ll2->SetTextSize(18);
  ll2->Draw();
  cname+=".pdf";
  c->Print(cname.data(),"pdf");
  string c2name = "c2_"+epname;
  TCanvas * c2 = new TCanvas(c2name.data(),c2name.data(),1000,600);
  c2->cd();
  rotOrig->SetMarkerColor(kRed);
  rotOrig->SetLineColor(kRed);
  rotOrig->SetMinimum(-15);
  rotOrig->SetMaximum(15);
  rotOrig->SetXTitle("Run #");
  rotOrig->SetYTitle("Angle Fluctuation (degrees)");
  rotOrig->SetMarkerStyle(24);
  rotOrig->SetMarkerSize(msize);
  rotOrig->GetXaxis()->SetLabelFont(43);
  rotOrig->GetXaxis()->SetLabelSize(22);
  rotOrig->GetYaxis()->SetLabelFont(43);
  rotOrig->GetYaxis()->SetLabelSize(22);
  rotOrig->GetYaxis()->SetTitleOffset(0.6);
  gPad->SetRightMargin(0.08);
  gPad->SetLeftMargin(0.075);
  rotOrig->Draw();
  rotFlat->SetMarkerColor(kBlue);
  rotFlat->SetLineColor(kBlue);
  //rotFlat->Draw("same");
  TLatex * detl = new TLatex(326500,12,epname.data());
  detl->SetTextFont(43);
  detl->SetTextSize(32);
  detl->Draw();
  for(int i = 0; i<51; i++) {
    l[i]=new TLine(iovs[i],15,iovs[i],17);
    l[i]->SetLineStyle(1);
    l[i]->Draw();
  }
  TLatex * ll3 = new TLatex(326400,15.5,"iovs:");
  ll3->SetTextFont(43);
  ll3->SetTextSize(18);
  ll3->Draw();
  string outname = c2name+".pdf";
  c2->Print(outname.data(),"pdf");
}
