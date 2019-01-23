TFile * f = new TFile("results.root","read");
int iovs[51]={326499, 326500, 326501, 326502, 326504, 326518, 326521, 326545, 
	    326547, 326550, 326574, 326575, 326576, 326580, 326585, 326586,
	    326587, 326617, 326620, 326627, 326635, 326790, 326822, 326854, 
	    326855, 326856, 326859, 326883, 326897, 326941, 326961, 326965, 
	    326998, 327004, 327078, 327123, 327125, 327147, 327148, 327173, 
	    327174, 327211, 327219, 327220, 327230, 327264, 327327, 327400, 
	    327403, 327424, 328000};

void plot(string epname = "HFm2"){
  string sqxOrig = "qxOrig_"+epname;
  string sqyOrig = "qyOrig_"+epname;
  string sqxFlat = "qxFlat_"+epname;
  string sqyFlat = "qyFlat_"+epname;
  TH1D * qxOrig = (TH1D *) f->Get(sqxOrig.data());
  TH1D * qyOrig = (TH1D *) f->Get(sqyOrig.data());
  TH1D * qxFlat = (TH1D *) f->Get(sqxFlat.data());
  TH1D * qyFlat = (TH1D *) f->Get(sqyFlat.data());
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
  qxOrig->SetMinimum(-0.2);
  qxOrig->SetMaximum(0.2);
  qxOrig->SetXTitle("Run #");
  qxOrig->SetYTitle("Offset");
  qxOrig->GetYaxis()->SetTitleOffset(0.8);
  qxOrig->Draw();
  qyOrig->Draw("same");
  TLatex * lep = new TLatex(326500,0.15,epname.data());
  lep->SetTextFont(43);
  lep->SetTextSize(32);
  lep->Draw();
  TLine * l[51];
  for(int i = 0; i<51; i++) {
    l[i]=new TLine(iovs[i],0.2,iovs[i],+0.25);
    l[i]->SetLineStyle(1);
    l[i]->Draw();
  }
  TLatex * ll = new TLatex(326400,0.205,"iovs:");
  ll->SetTextFont(43);
  ll->SetTextSize(18);
  ll->Draw();
  c->cd(2);
  gPad->SetLeftMargin(0.1);
  qxFlat->SetMinimum(-0.2);
  qxFlat->SetMaximum(0.2);
  qxFlat->SetXTitle("Run #");
  qxFlat->SetYTitle("Offset");
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
}
