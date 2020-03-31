bool Fill_N(int anal,int n, int bin, TH2F *qxtrk_, TH2F * qytrk_, TH2F * qcnt_, double Ax, double Ay, double Bx, double By, double Cx, double Cy, double wA, double wB, double wC, double thA, double thB, double thC){
  bool stat = false;
  if(pow(Ax,2)+pow(Ay,2) < 1e-6) return stat;
  if(pow(Bx,2)+pow(By,2) < 1e-6) return stat;
  if(pow(Cx,2)+pow(Cy,2) < 1e-6) return stat;
  for(int i = 1; i<=qxtrk_->GetNbinsX(); i++) {
    for(int j = 1; j<=qxtrk_->GetNbinsY(); j++) {
      if(qxtrk_->GetBinContent(i,j)!=0) {
	double th = TMath::ATan2(qytrk_->GetBinContent(i,j),qxtrk_->GetBinContent(i,j));
	qanal[anal].rA[bin][0]->SetBinContent(i,j,qanal[anal].rA[bin][0]->GetBinContent(i,j)+TMath::Cos(n*(th-thA)));
	qanal[anal].rB[bin][0]->SetBinContent(i,j,qanal[anal].rB[bin][0]->GetBinContent(i,j)+TMath::Cos(n*(th-thB)));
	qanal[anal].rAcnt[bin][0]->SetBinContent(i,j,qanal[anal].rAcnt[bin][0]->GetBinContent(i,j)+1);
	qanal[anal].rBcnt[bin][0]->SetBinContent(i,j,qanal[anal].rBcnt[bin][0]->GetBinContent(i,j)+1);
      }
    }
  }
  qanal[anal].qA[bin][0]->Add(qxtrk_,Ax);
  qanal[anal].qA[bin][0]->Add(qytrk_,Ay);
  qanal[anal].qB[bin][0]->Add(qxtrk_,Bx);
  qanal[anal].qB[bin][0]->Add(qytrk_,By);
  qanal[anal].wnA[bin][0]->Add(qcnt_,wA);
  qanal[anal].wnB[bin][0]->Add(qcnt_,wB);
  qanal[anal].qBA[bin][0]->Fill(0.,Bx*Ax + By*Ay);
  qanal[anal].qCA[bin][0]->Fill(0.,Cx*Ax + Cy*Ay);
  qanal[anal].qCB[bin][0]->Fill(0.,Cx*Bx + Cy*By);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,wB*wA);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,wC*wA);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,wC*wB);
  qanal[anal].wA[bin][0]->Fill(wA);
  qanal[anal].wB[bin][0]->Fill(wB);
  qanal[anal].wC[bin][0]->Fill(wC);
  qanal[anal].rBA[bin][0]->Fill(0.,TMath::Cos(n*(thB-thA)));
  qanal[anal].rCA[bin][0]->Fill(0.,TMath::Cos(n*(thC-thA)));
  qanal[anal].rCB[bin][0]->Fill(0.,TMath::Cos(n*(thC-thB)));
  qanal[anal].rBAcnt[bin][0]->Fill(0.);
  qanal[anal].rCAcnt[bin][0]->Fill(0.);
  qanal[anal].rCBcnt[bin][0]->Fill(0.);
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,Ax);
  qanal[anal].qA[bin][j]->Add(qytrk_,Ay);
  qanal[anal].qB[bin][j]->Add(qxtrk_,Bx);
  qanal[anal].qB[bin][j]->Add(qytrk_,By);
  qanal[anal].wnA[bin][j]->Add(qcnt_,wA);
  qanal[anal].wnB[bin][j]->Add(qcnt_,wB);
  qanal[anal].qBA[bin][j]->Fill(0.,Bx*Ax + By*Ay);
  qanal[anal].qCA[bin][j]->Fill(0.,Cx*Ax + Cy*Ay);
  qanal[anal].qCB[bin][j]->Fill(0.,Cx*Bx + Cy*By);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,wB*wA);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,wC*wA);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,wC*wB);
  qanal[anal].wA[bin][j]->Fill(wA);
  qanal[anal].wB[bin][j]->Fill(wB);
  qanal[anal].wC[bin][j]->Fill(wC);
  qanal[anal].rBA[bin][j]->Fill(0.,TMath::Cos(n*(thB-thA)));
  qanal[anal].rCA[bin][j]->Fill(0.,TMath::Cos(n*(thC-thA)));
  qanal[anal].rCB[bin][j]->Fill(0.,TMath::Cos(n*(thC-thB)));
  qanal[anal].rBAcnt[bin][j]->Fill(0.);
  qanal[anal].rCAcnt[bin][j]->Fill(0.);
  qanal[anal].rCBcnt[bin][j]->Fill(0.);
  stat = true;
  for(int i = 1; i<=qxtrk_->GetNbinsX(); i++) {
    for(int k = 1; k<=qxtrk_->GetNbinsY(); k++) {
      if(qxtrk_->GetBinContent(i,k)!=0) {
	double th = TMath::ATan2(qytrk_->GetBinContent(i,k),qxtrk_->GetBinContent(i,k));
	qanal[anal].rA[bin][j]->SetBinContent(i,k,qanal[anal].rA[bin][j]->GetBinContent(i,k)+TMath::Cos(n*(th-thA)));
	qanal[anal].rB[bin][j]->SetBinContent(i,k,qanal[anal].rB[bin][j]->GetBinContent(i,k)+TMath::Cos(n*(th-thB)));
	qanal[anal].rAcnt[bin][j]->SetBinContent(i,k,qanal[anal].rAcnt[bin][0]->GetBinContent(i,k)+1);
	qanal[anal].rBcnt[bin][j]->SetBinContent(i,k,qanal[anal].rBcnt[bin][0]->GetBinContent(i,k)+1);
      }
    }
  }
  return stat;
}
