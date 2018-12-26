#include "../EPCalib/HiEvtPlaneList.h"
using namespace hi;
TFile * fin;
int checkfiles[4]={HFm2,HFp2,trackmid2,Castor2};
TH1D * hflat[4][50];
static const int minrun = 326381;
static const int maxrun = 327564;
static const int nruns = 1184;
TH1D * runs = new TH1D("runs","runs",nruns,minrun,maxrun);
int iruns[nruns];
int goodruns[50];
int ngoodruns = 0;
TTree * tr;
void CheckFlattening(string fname="../save/EP_PbPb2018_1_326886.root") {
  for(int i = 0; i<4; i++) {
    for(int j = 0; j< 50; j++) {
      //      string rawname = "raw_"+checkfiles[i]+"_"+to_string(j);
      string flatname = "flat_"+EPNames[checkfiles[i]]+"_"+to_string(j);
      //      TH1D * hraw[i][j] = new TH1D(rawname.data(),rawname.data(),200,-3,3);
      hflat[i][j] = new TH1D(flatname.data(),flatname.data(),100,-2,2);

    }
  }
  float EPcent;
  float EPvtx;
  unsigned int EPntrk;
  unsigned int EPrun;
  float EPAngs[NumEPNames];
  int ptcnt[NumEPNames];

  fin = new TFile(fname.data());
  if(fin->IsZombie()) {
    cout<<fname<<"  does not exist or is corrupted."<<endl;
    return 0;
  }
  tr = (TTree *) fin->Get("EPtree");
  tr->SetBranchAddress("Cent",    &EPcent);
  tr->SetBranchAddress("Vtx",     &EPvtx);
  tr->SetBranchAddress("Run",     &EPrun);
  tr->SetBranchAddress("EPAngs", &EPAngs);
  
  int nentries = tr->GetEntries();
  cout<<nentries<<endl;
  memset(iruns,0,sizeof(iruns));
  for(int i = 0; i< nentries/50; i++) {
    tr->GetEntry(i);
    ++iruns[EPrun-minrun];
    runs->Fill(EPrun);
  }
  for(int i = minrun; i<=maxrun; i++) {
    if(iruns[i-minrun]>0) goodruns[ngoodruns++]=i;
  }

  for(int i = 0; i< nentries; i++) {
    tr->GetEntry(i);
    bool sel = false;
    int nsel = -1;
    for(int j = 0; j<ngoodruns; j++) {
      if(EPrun==goodruns[j]) {
	sel = true;
	nsel = j;
	break;
      }
    }    
    if(sel) {
      for(int k=0; k<4; k++) {
	hflat[k][nsel]->Fill(EPAngs[checkfiles[k]]);
      }
    }
  }

  TCanvas * can[50];
  for(int i = 0; i< ngoodruns; i++) {
    string canname = "";
    if(fname.find("general")!=std::string::npos) {
      canname+="general_";
    } else if(fname.find("pixel")!=std::string::npos) {
      canname+="pixel_";
    } else canname+="unknown_";
    canname += to_string(goodruns[i]);
   
    can[i] = new TCanvas(canname.data(),canname.data(),1100,1000);
    can[i]->Divide(2,2);
    for(int j = 0; j< 4; j++) {
      can[i]->cd(j+1);
      hflat[j][i]->Draw();
      TLatex * name = new TLatex(-0.5,0.5*hflat[j][i]->GetMaximum(),EPNames[checkfiles[j]].data());
      name->Draw();
      string fname = canname+"fit"+to_string(j);
      TF1 * f = new TF1(fname.data(),"pol0",-1.5,1.5);
      hflat[j][i]->Fit(f,"qrn");
      cout<<goodruns[i]<<"\t"<<EPNames[checkfiles[j]]<<"\t"<<f->GetChisquare()/f->GetNDF()<<endl;
      TLatex * rchi = new TLatex(-0.5,0.4*hflat[j][i]->GetMaximum(),Form("#chi^{2}/NDF = %5.2f",f->GetChisquare()/f->GetNDF()));
      rchi->Draw();
      float rchisqr = f->GetChisquare()/f->GetNDF();
      std::ostringstream out;
      out.precision(2);
      out<<std::fixed<<rchisqr;
      canname+="_"+out.str();
    }
    canname= "checkplots/"+canname+".pdf";
    can[i]->Print(canname.data(),"pdf");
  }
  fin->Close();
  
}
