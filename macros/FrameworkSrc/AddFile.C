#include <ctime>

bool Framework::AddFile(){
  //clock_t start = clock();
  char buf[120];
  if(fgets(buf,120,flist)==NULL) return false;
  buf[strlen(buf)-1]=0;
  string infile = buf;
  //cout<<infile.data()<<endl;
  std::unique_ptr<TFile> tf(new TFile(infile.data(),"read"));
  if(tf->IsZombie())                 {
   cout<<"ZOMBIE:    " <<infile.data()<<endl; 
   string remove = "rm "+infile;
   system(remove.data());
   return true;
  }

  int minCent = r[0].minCent;
  int maxCent = r[0].maxCent;
  int minCentBin = hcentbins->FindBin(minCent);
  int maxCentBin = hcentbins->FindBin(maxCent-0.1);
  cout<<"centrality range: "<<minCent<<" - "<<maxCent<<"\t  bins: "<<minCentBin<<"\t"<<maxCentBin<<endl;
  TH2F * ptspec=0;
  TH2F * ptspecCnt=0;
  TH1D * ntrk=0;
  TH2F * qA[11]={0};
  TH2F * qB[11]={0};
  TH2F * wnA[11]={0};
  TH2F * wnB[11]={0};
  TH1D * qBA[11]={0};
  TH1D * qCA[11]={0};
  TH1D * qCB[11]={0};
  TH1D * qBAcnt[11]={0};
  TH1D * qCAcnt[11]={0};
  TH1D * qCBcnt[11]={0};
  for(int icent = minCentBin; icent<= maxCentBin; icent++) {
    int cmin = hcentbins->GetBinLowEdge(icent);
    int cmax = hcentbins->GetBinLowEdge(icent)+hcentbins->GetBinWidth(icent);
    if(icent==minCentBin) {
      ptspec = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/ptav",cmin,cmax));
      ptspec->SetDirectory(0);
      ptspecCnt = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/ptcnt",cmin,cmax));
      ptspecCnt->SetDirectory(0);
      qA[0] = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qA",cmin,cmax,r[0].order));
      qA[0]->SetDirectory(0);
      qB[0] = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qB",cmin,cmax,r[0].order));
      qB[0]->SetDirectory(0);
      wnA[0] = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/wnA",cmin,cmax,r[0].order));
      wnA[0]->SetDirectory(0);
      wnB[0] = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/wnB",cmin,cmax,r[0].order));
      wnB[0]->SetDirectory(0);
      qBA[0] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qBA",cmin,cmax,r[0].order));
      qBA[0]->SetDirectory(0);
      qCA[0] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qCA",cmin,cmax,r[0].order));
      qCA[0]->SetDirectory(0);
      qCB[0] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qCB",cmin,cmax,r[0].order));
      qCB[0]->SetDirectory(0);
      qBAcnt[0] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qBAcnt",cmin,cmax,r[0].order));
      qBAcnt[0]->SetDirectory(0);
      qCAcnt[0] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qCAcnt",cmin,cmax,r[0].order));
      qCAcnt[0]->SetDirectory(0);
      qCBcnt[0] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qCBcnt",cmin,cmax,r[0].order));
      qCBcnt[0]->SetDirectory(0);
      for(int isub = 1; isub<=10; isub++) {
	qA[isub] = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qA_%d",cmin,cmax,r[0].order,isub));
	qA[isub]->SetDirectory(0);
	qB[isub] = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qB_%d",cmin,cmax,r[0].order,isub));
	qB[isub]->SetDirectory(0);
	wnA[isub] = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/wnA_%d",cmin,cmax,r[0].order,isub));
	wnA[isub]->SetDirectory(0);
	wnB[isub] = (TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/wnB_%d",cmin,cmax,r[0].order,isub));
	wnB[isub]->SetDirectory(0);
	qBA[isub] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qBA_%d",cmin,cmax,r[0].order,isub));
	qBA[isub]->SetDirectory(0);
	qCA[isub] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qCA_%d",cmin,cmax,r[0].order,isub));
	qCA[isub]->SetDirectory(0);
	qCB[isub] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qCB_%d",cmin,cmax,r[0].order,isub));
	qCB[isub]->SetDirectory(0);
	qBAcnt[isub] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qBAcnt_%d",cmin,cmax,r[0].order,isub));
	qBAcnt[isub]->SetDirectory(0);
	qCAcnt[isub] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qCAcnt_%d",cmin,cmax,r[0].order,isub));
	qCAcnt[isub]->SetDirectory(0);
	qCBcnt[isub] = (TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qCBcnt_%d",cmin,cmax,r[0].order,isub));
	qCBcnt[isub]->SetDirectory(0);
      }
    } else {
      ptspec->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/ptav",cmin,cmax)));
      ptspecCnt->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/ptcnt",cmin,cmax)));
      qA[0]->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qA",cmin,cmax,r[0].order)));
      qB[0]->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qB",cmin,cmax,r[0].order)));
      wnA[0]->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/wnA",cmin,cmax,r[0].order)));
      wnB[0]->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/wnB",cmin,cmax,r[0].order)));
      qBA[0]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qBA",cmin,cmax,r[0].order)));
      qCA[0]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qCA",cmin,cmax,r[0].order)));
      qCB[0]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qCB",cmin,cmax,r[0].order)));
      qBAcnt[0]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qBAcnt",cmin,cmax,r[0].order)));
      qCAcnt[0]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qCAcnt",cmin,cmax,r[0].order)));
      qCBcnt[0]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/qCBcnt",cmin,cmax,r[0].order)));
      for(int isub = 1; isub<=10; isub++) {
	  qA[isub]->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qA_%d",cmin,cmax,r[0].order,isub)));
	  qB[isub]->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qB_%d",cmin,cmax,r[0].order,isub)));
	  wnA[isub]->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/wnA_%d",cmin,cmax,r[0].order,isub)));
	  wnB[isub]->Add((TH2F *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/wnB_%d",cmin,cmax,r[0].order,isub)));
	  qBA[isub]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qBA_%d",cmin,cmax,r[0].order,isub)));
	  qCA[isub]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qCA_%d",cmin,cmax,r[0].order,isub)));
	  qCB[isub]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qCB_%d",cmin,cmax,r[0].order,isub)));
	  qBAcnt[isub]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qBAcnt_%d",cmin,cmax,r[0].order,isub)));
	  qCAcnt[isub]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qCAcnt_%d",cmin,cmax,r[0].order,isub)));
	  qCBcnt[isub]->Add((TH1D *) tf->Get(Form("vnanalyzer/Harmonics/%d_%d/N%d/SubEvents/qCBcnt_%d",cmin,cmax,r[0].order,isub)));
      }
      hcent->Add((TH1D *) tf->Get("vnanalyzer/cent"));
    }
  }
  double etamin1 = r[0].minEta[0];
  double etamax1 = r[0].maxEta[0];
  int ietamin1 = ptspecCnt->GetYaxis()->FindBin(etamin1+0.01);
  int ietamax1 = ptspecCnt->GetYaxis()->FindBin(etamax1-0.01);
  double etamin2 = r[0].minEta[1];
  double etamax2 = r[0].maxEta[1];
    int ietamin2 = ptspecCnt->GetYaxis()->FindBin(etamin2+0.01);
    int ietamax2 = ptspecCnt->GetYaxis()->FindBin(etamax2-0.01);
    if(r[0].spec_NegEta == 0) {
      r[0].spec_NegEta = (TH1D *) ptspecCnt->ProjectionX(Form("spec1d_%d_%d_NegEta",minCent,maxCent),ietamin1,ietamax1);
      r[0].spec_NegEta->SetDirectory(0);
      r[0].spec_NegEta->SetMarkerStyle(20);
      r[0].spec_NegEta->SetMarkerColor(kBlue);
    } else {
      r[0].spec_NegEta->Add((TH1D *) ptspecCnt->ProjectionX(Form("spec1d_%d_%d_NegEta",minCent,maxCent),ietamin1,ietamax1));
    }
    if(r[0].spec_PosEta == 0 ) {
      r[0].spec_PosEta = (TH1D *) ptspecCnt->ProjectionX(Form("spec1d_%d_%d_PosEta",minCent,maxCent),ietamin2,ietamax2);
      r[0].spec_PosEta->SetDirectory(0);
      r[0].spec_PosEta->SetMarkerStyle(25);
      r[0].spec_PosEta->SetMarkerColor(kRed);
    } else {
      r[0].spec_PosEta->Add((TH1D *) ptspecCnt->ProjectionX(Form("spec1d_%d_%d_PosEta",minCent,maxCent),ietamin2,ietamax2));
    }
    tf->Close();
  
  for(int i = 0; i<nrange; i++) {
    for(int j = r[i].minEtaBin[0]; j<= r[i].maxEtaBin[0]; j++) {
      for(int k = r[i].minPtBin; k<=r[i].maxPtBin; k++) {
  	r[i].qn[0]+=qA[0]->GetBinContent(k,j);
  	r[i].qne[0]+=pow(qA[0]->GetBinError(k,j),2);
  	r[i].wnA[0] += wnA[0]->GetBinContent(k,j);
  	r[i].qAB[0] += qBA[0]->GetBinContent(1);
  	r[i].qAC[0] += qCA[0]->GetBinContent(1);
  	r[i].qBC[0] += qCB[0]->GetBinContent(1);
  	r[i].wAB[0] += qBAcnt[0]->GetBinContent(1);
  	r[i].wAC[0] += qCAcnt[0]->GetBinContent(1);
  	r[i].wBC[0] += qCBcnt[0]->GetBinContent(1);
	
  	r[i].pt[0] += ptspec->GetBinContent(k,j);
  	r[i].ptcnt[0] += ptspecCnt->GetBinContent(k,j);
	
  	for(int m = 1; m<11; m++) {  
  	  r[i].qnSub[0][m-1]+=qA[m]->GetBinContent(k,j);
  	  r[i].wnASub[0][m-1] += wnA[m]->GetBinContent(k,j);
  	}
      }
    }
    //////////////////
    for(int j = r[i].minEtaBin[1]; j<= r[i].maxEtaBin[1]; j++) {
      for(int k = r[i].minPtBin; k<=r[i].maxPtBin; k++) {
  	r[i].qn[1]+=qB[0]->GetBinContent(k,j);
  	r[i].qne[1]+=pow(qB[0]->GetBinError(k,j),2);
  	r[i].wnA[1] += wnB[0]->GetBinContent(k,j);
  	r[i].qAB[1] += qBA[0]->GetBinContent(1);
  	r[i].qAC[1] += qCB[0]->GetBinContent(1);
  	r[i].qBC[1] += qCA[0]->GetBinContent(1);
  	r[i].wAB[1] += qBAcnt[0]->GetBinContent(1);
  	r[i].wAC[1] += qCBcnt[0]->GetBinContent(1);
  	r[i].wBC[1] += qCAcnt[0]->GetBinContent(1);
	
  	r[i].pt[1] += ptspec->GetBinContent(k,j);
  	r[i].ptcnt[1] += ptspecCnt->GetBinContent(k,j);
	
  	for(int m = 1; m<11; m++) {  
  	  r[i].qnSub[1][m-1]+=qB[m]->GetBinContent(k,j);
  	  r[i].wnASub[1][m-1] += wnB[m]->GetBinContent(k,j);
  	}
      }
    }
  }
  if(ptspec) ptspec->Delete("*/*");
  if(ptspecCnt) ptspecCnt->Delete("*/*");
  if(ntrk) ntrk->Delete("*/*");
  for(int k = 0; k<11; k++) {
    if(qA[k]) {
      qA[k]->Delete("*/*");
      qA[k] = 0;
    }
    if(qB[k]) {
      qB[k]->Delete("*/*");
      qB[k] = 0;
    }
    if(wnA[k]) {
      wnA[k]->Delete("*/*");
      wnA[k] = 0;
    }
    if(wnB[k]) {
      wnB[k]->Delete("*/*");
      wnB[k] = 0;
    }
    if(qBA[k]) {
      qBA[k]->Delete("*/*");
      qBA[k] = 0;
    }
    if(qCA[k]) {
      qCA[k]->Delete("*/*");
      qCA[k] = 0;
    }
    if(qCB[k]) {
      qCB[k]->Delete("*/*");
      qCB[k] = 0;
    }
    if(qBAcnt[k]) {
      qBAcnt[k]->Delete("*/*");
      qBAcnt[k] = 0;
    }
    if(qCAcnt[k]) {
      qCAcnt[k]->Delete("*/*");
      qCAcnt[k]=0;
    }
    if(qCBcnt[k]) {
      qCBcnt[k]->Delete("*/*");
      qCBcnt[k] = 0;
    }
  }
  //clock_t end = clock();
  //double elapsed_secs = double(end-start)/CLOCKS_PER_SEC;
  //cout<<"elapsed time with file (s) = "<<elapsed_secs<<endl;
  
  return true;
}
