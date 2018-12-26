



void Fill_N(int anal,double order, int bin, TH2D *qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double Ax, double Ay, double Bx, double By, double Cx, double Cy, double wA, double wB, double wC){
  if(pow(Ax,2)+pow(Ay,2) < 1e-6) return;
  if(pow(Bx,2)+pow(By,2) < 1e-6) return;
  if(pow(Cx,2)+pow(Cy,2) < 1e-6) return;
  
  qanal[anal].thA[bin][0]->Fill(TMath::ATan2(Ay,Ax)/order);
  qanal[anal].thB[bin][0]->Fill(TMath::ATan2(By,Bx)/order);
  qanal[anal].thC[bin][0]->Fill(TMath::ATan2(Cy,Cx)/order);
  qanal[anal].thA2[bin][0]->Fill(0);
  qanal[anal].thB2[bin][0]->Fill(0);
  qanal[anal].thC2[bin][0]->Fill(0);
  qanal[anal].thn[bin][0]->Fill(TMath::ATan2(qytrk_->GetBinContent(5,9),qxtrk_->GetBinContent(5,9))/order);


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
}

void Fill_N112(int anal, int bin, int a1, int b1, int c1, int a2, int b2, int c2, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a1x = qx[a1];
  double a1y = qy[a1];
  double a1w = sumw[a1];
  double b1x = qx[b1];
  double b1y = qy[b1];
  double b1w = sumw[b1];
  double c1x = qx[c1];
  double c1y = qy[c1];
  double c1w = sumw[c1];

  double a2x = qx[a2];
  double a2y = qy[a2];
  double a2w = sumw[a2];
  double b2x = qx[b2];
  double b2y = qy[b2];
  double b2w = sumw[b2];
  double c2x = qx[c2];
  double c2y = qy[c2];
  double c2w = sumw[c2];


  if(pow(a1x,2)+pow(a1y,2) < 1e-12) return;
  if(pow(b1x,2)+pow(b1y,2) < 1e-12) return;
  if(pow(c1x,2)+pow(c1y,2) < 1e-12) return;
  if(pow(a2x,2)+pow(a2y,2) < 1e-12) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-12) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-12) return;
  
  qanal[anal].thA[bin][0]->Fill(TMath::ATan2(a1y,a1x)/EPOrder[a1]);
  qanal[anal].thB[bin][0]->Fill(TMath::ATan2(b1y,b1x)/EPOrder[b1]);
  qanal[anal].thC[bin][0]->Fill(TMath::ATan2(c1y,c1x)/EPOrder[c1]);
  qanal[anal].thA2[bin][0]->Fill(TMath::ATan2(a2y,a2x)/EPOrder[a2]);
  qanal[anal].thB2[bin][0]->Fill(TMath::ATan2(b2y,b2x)/EPOrder[b2]);
  qanal[anal].thC2[bin][0]->Fill(TMath::ATan2(c2y,c2x)/EPOrder[c2]);
  qanal[anal].thn[bin][0]->Fill(TMath::ATan2(qytrk_->GetBinContent(5,9),qxtrk_->GetBinContent(5,9))/EPOrder[a1]);

  qanal[anal].qA[bin][0]->Add(qxtrk_,  a1x*a2x+a1y*a2y);
  qanal[anal].qA[bin][0]->Add(qytrk_, -a1y*a2x+a1x*a2y);
  qanal[anal].qB[bin][0]->Add(qxtrk_,  b1x*a2x+b1y*a2y);
  qanal[anal].qB[bin][0]->Add(qytrk_, -b1y*a2x+b1x*a2y);
  qanal[anal].wnA[bin][0]->Add(qcnt_,a1w*a2w);
  qanal[anal].wnB[bin][0]->Add(qcnt_,b1w*a2w);
  qanal[anal].qBA[bin][0]->Fill(0., a1x*b1x + a1y*b1y);
  qanal[anal].qCA[bin][0]->Fill(0., a1x*c1x + a1y*c1y); 
  qanal[anal].qCB[bin][0]->Fill(0., b1x*c1x + b1y*c1y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b1w*a1w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,a1w*c1w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,b1w*c1w);
  qanal[anal].qBA2[bin][0]->Fill(0., a2x*b2x + a2y*b2y);
  qanal[anal].qCA2[bin][0]->Fill(0., a2x*c2x + a2y*c2y); 
  qanal[anal].qCB2[bin][0]->Fill(0., b2x*c2x + b2y*c2y);
  qanal[anal].qBAcnt2[bin][0]->Fill(0.,b2w*a2w);
  qanal[anal].qCAcnt2[bin][0]->Fill(0.,a2w*c2w);
  qanal[anal].qCBcnt2[bin][0]->Fill(0.,b2w*c2w);
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,  a1x*a2x+a1y*a2y);
  qanal[anal].qA[bin][j]->Add(qytrk_, -a1y*a2x+a1x*a2y);
  qanal[anal].qB[bin][j]->Add(qxtrk_,  b1x*a2x+b1y*a2y);
  qanal[anal].qB[bin][j]->Add(qytrk_, -b1y*a2x+b1x*a2y);
  qanal[anal].wnA[bin][j]->Add(qcnt_,a1w*a2w);
  qanal[anal].wnB[bin][j]->Add(qcnt_,b1w*a2w);
  qanal[anal].qBA[bin][j]->Fill(0., a1x*b1x + a1y*b1y);
  qanal[anal].qCA[bin][j]->Fill(0., a1x*c1x + a1y*c1y); 
  qanal[anal].qCB[bin][j]->Fill(0., b1x*c1x + b1y*c1y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b1w*a1w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,a1w*c1w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,b1w*c1w);
  qanal[anal].qBA2[bin][j]->Fill(0., a2x*b2x + a2y*b2y);
  qanal[anal].qCA2[bin][j]->Fill(0., a2x*c2x + a2y*c2y); 
  qanal[anal].qCB2[bin][j]->Fill(0., b2x*c2x + b2y*c2y);
  qanal[anal].qBAcnt2[bin][j]->Fill(0.,b2w*a2w);
  qanal[anal].qCAcnt2[bin][j]->Fill(0.,a2w*c2w);
  qanal[anal].qCBcnt2[bin][j]->Fill(0.,b2w*c2w);
}

void Fill_N112Combined(int anal, int bin, int a1, int b1, int c1, int a2, int b2, int c2,TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){

  double a1x = qx[a1];
  double a1y = qy[a1];
  double b1x = qx[b1];
  double b1y = qy[b1];
  double c1x = qx[c1];
  double c1y = qy[c1];
  
  double a2x = qx[a2];
  double a2y = qy[a2];
  double b2x = qx[b2];
  double b2y = qy[b2];
  double c2x = qx[c2];
  double c2y = qy[c2];
  
  if(pow(a1x,2)+pow(a1y,2) < 1e-8) return;
  if(pow(b1x,2)+pow(b1y,2) < 1e-8) return;
  if(pow(c1x,2)+pow(c1y,2) < 1e-8) return;
  if(pow(a2x,2)+pow(a2y,2) < 1e-8) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-8) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-8) return;
  
  a2x/=sqrt(pow(a2x,2)+pow(a2y,2));
  a2y/=sqrt(pow(a2x,2)+pow(a2y,2));
  b2x/=sqrt(pow(b2x,2)+pow(b2y,2));
  b2y/=sqrt(pow(b2x,2)+pow(b2y,2));
  c2x/=sqrt(pow(c2x,2)+pow(c2y,2));
  c2y/=sqrt(pow(c2x,2)+pow(c2y,2));
  
  double psi1A = TMath::ATan2(a1y,a1x);
  double psi1B = TMath::ATan2(b1y,b1x)+TMath::Pi();
  double psi2A = TMath::ATan2(a2y,a2x)/2.;
  double psi2B = TMath::ATan2(b2y,b2x)/2.;
  double psi2C = TMath::ATan2(c2y,c2x)/2.;
  qanal[anal].qA[bin][0]->Add(qxtrk_, TMath::Cos(psi1A-2*psi2A));
  qanal[anal].qA[bin][0]->Add(qytrk_, -TMath::Sin(psi1A-2*psi2A));
  qanal[anal].qB[bin][0]->Add(qxtrk_, TMath::Cos(psi1B-2*psi2A));
  qanal[anal].qB[bin][0]->Add(qytrk_, -TMath::Sin(psi1B-2*psi2A));
  qanal[anal].wnA[bin][0]->Add(qcnt_);
  qanal[anal].wnB[bin][0]->Add(qcnt_);
  qanal[anal].qBA[bin][0]->Fill(0., TMath::Cos(psi1A+psi1B-2*psi2A));
  qanal[anal].qCA[bin][0]->Fill(0.); 
  qanal[anal].qCB[bin][0]->Fill(0.);
  qanal[anal].qBAcnt[bin][0]->Fill(0.);
  qanal[anal].qCAcnt[bin][0]->Fill(0.);
  qanal[anal].qCBcnt[bin][0]->Fill(0.);
  qanal[anal].qBA2[bin][0]->Fill(0., TMath::Cos(psi2A-psi2B));
  qanal[anal].qCA2[bin][0]->Fill(0., TMath::Cos(psi2A-psi2C)); 
  qanal[anal].qCB2[bin][0]->Fill(0., TMath::Cos(psi2B-psi2C));
  qanal[anal].qBAcnt2[bin][0]->Fill(0.);
  qanal[anal].qCAcnt2[bin][0]->Fill(0.);
  qanal[anal].qCBcnt2[bin][0]->Fill(0.);
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_, TMath::Cos(psi1A-2*psi2A));
  qanal[anal].qA[bin][j]->Add(qytrk_, -TMath::Sin(psi1A-2*psi2A));
  qanal[anal].qB[bin][j]->Add(qxtrk_, TMath::Cos(psi1B-2*psi2A));
  qanal[anal].qB[bin][j]->Add(qytrk_, -TMath::Sin(psi1B-2*psi2A));
  qanal[anal].wnA[bin][j]->Add(qcnt_);
  qanal[anal].wnB[bin][j]->Add(qcnt_);
  qanal[anal].qBA[bin][j]->Fill(0., TMath::Cos(psi1A+psi1B-2*psi2A));
  qanal[anal].qCA[bin][j]->Fill(0.); 
  qanal[anal].qCB[bin][j]->Fill(0.);
  qanal[anal].qBAcnt[bin][j]->Fill(0.);
  qanal[anal].qCAcnt[bin][j]->Fill(0.);
  qanal[anal].qCBcnt[bin][j]->Fill(0.);
  qanal[anal].qBA2[bin][j]->Fill(0., TMath::Cos(psi2A-psi2B));
  qanal[anal].qCA2[bin][j]->Fill(0., TMath::Cos(psi2A-psi2C)); 
  qanal[anal].qCB2[bin][j]->Fill(0., TMath::Cos(psi2B-psi2C));
  qanal[anal].qBAcnt2[bin][j]->Fill(0.);
  qanal[anal].qCAcnt2[bin][j]->Fill(0.);
  qanal[anal].qCBcnt2[bin][j]->Fill(0.);
}

void Fill_N42(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double Ax = qx[HFp2];
  double Ay = qy[HFp2];
  double wA = sumw[HFp2];

  double Bx = qx[HFm2];
  double By = qy[HFm2];
  double wB = sumw[HFm2];

  double Cx = qx[trackmid2];
  double Cy = qy[trackmid2];
  double wC = sumw[trackmid2];

  if(pow(Ax,2)+pow(Ay,2) < 1e-6) return;
  if(pow(Bx,2)+pow(By,2) < 1e-6) return;
  if(pow(Cx,2)+pow(Cy,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,pow(Ax,2)-pow(Ay,2));
  qanal[anal].qA[bin][0]->Add(qytrk_,2*Ax*Ay);
  qanal[anal].qB[bin][0]->Add(qxtrk_,pow(Bx,2)-pow(By,2));
  qanal[anal].qB[bin][0]->Add(qytrk_,2*Bx*By);
  qanal[anal].wnA[bin][0]->Add(qcnt_,pow(wA,2));
  qanal[anal].wnB[bin][0]->Add(qcnt_,pow(wB,2));
  qanal[anal].qBA[bin][0]->Fill(0.,pow(Ax,2)*(pow(Bx,2)-pow(By,2)) + pow(Ay,2)*(pow(By,2)-pow(Bx,2))+4*Ax*Ay*Bx*By);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(Ax,2)*(pow(Cx,2)-pow(Cy,2)) + pow(Ay,2)*(pow(Cy,2)-pow(Cx,2))+4*Ax*Ay*Cx*Cy);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(Bx,2)*(pow(Cx,2)-pow(Cy,2)) + pow(By,2)*(pow(Cy,2)-pow(Cx,2))+4*Bx*By*Cx*Cy);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wA*wB,2));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wA*wC,2));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wB*wC,2));
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,pow(Ax,2)-pow(Ay,2));
  qanal[anal].qA[bin][j]->Add(qytrk_,2*Ax*Ay);
  qanal[anal].qB[bin][j]->Add(qxtrk_,pow(Bx,2)-pow(By,2));
  qanal[anal].qB[bin][j]->Add(qytrk_,2*Bx*By);
  qanal[anal].wnA[bin][j]->Add(qcnt_,pow(wA,2));
  qanal[anal].wnB[bin][j]->Add(qcnt_,pow(wB,2));
  qanal[anal].qBA[bin][j]->Fill(0.,pow(Ax,2)*(pow(Bx,2)-pow(By,2)) + pow(Ay,2)*(pow(By,2)-pow(Bx,2))+4*Ax*Ay*Bx*By);
  qanal[anal].qCA[bin][j]->Fill(0.,pow(Ax,2)*(pow(Cx,2)-pow(Cy,2)) + pow(Ay,2)*(pow(Cy,2)-pow(Cx,2))+4*Ax*Ay*Cx*Cy);
  qanal[anal].qCB[bin][j]->Fill(0.,pow(Bx,2)*(pow(Cx,2)-pow(Cy,2)) + pow(By,2)*(pow(Cy,2)-pow(Cx,2))+4*Bx*By*Cx*Cy);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wA*wB,2));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wA*wC,2));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wB*wC,2));
}

void Fill_N42A(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp2];
  double ay = qy[HFp2];
  double aw = sumw[HFp2];

  double bx = qx[HFm2];
  double by = qy[HFm2];
  double bw = sumw[HFm2];

  double cx = qx[trackmid2];
  double cy = qy[trackmid2];
  double cw = sumw[trackmid2];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][0]->Add(qcnt_,aw*bw);
  qanal[anal].wnB[bin][0]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(aw*bw,2));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(aw*cw,2));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(bw*cw,2));

  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][j]->Add(qcnt_,aw*aw);
  qanal[anal].wnB[bin][j]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by);
  qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy);
  qanal[anal].qCB[bin][j]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(aw*bw,2));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(aw*cw,2));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(bw*cw,2));

}
void Fill_N42B(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp2];
  double ay = qy[HFp2];
  double aw = sumw[HFp2];

  double bx = qx[HFm2];
  double by = qy[HFm2];
  double bw = sumw[HFm2];

  double cx = qx[trackmid2];
  double cy = qy[trackmid2];
  double cw = sumw[trackmid2];

  double dx = qx[trackp210];
  double dy = qy[trackp210];
  double dw = sumw[trackp210];

  double ex = qx[trackm210];
  double ey = qy[trackm210];
  double ew = sumw[trackm210];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;
  if(pow(dx,2)+pow(dy,2) < 1e-6) return;
  if(pow(ex,2)+pow(ey,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][0]->Add(qcnt_,aw*bw);
  qanal[anal].wnB[bin][0]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][0]->Fill(0.,ax*bx*dx*ex - ay*by*dx*ex + ay*bx*dy*ex + ax*by*dy*ex + ay*bx*dx*ey + ax*by*dx*ey - ax*bx*dy*ey + ay*by*dy*ey);
  qanal[anal].qCA[bin][0]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2));
  qanal[anal].qCB[bin][0]->Fill(0.,pow(cx,2)*dx*ex - pow(cy,2)*dx*ex + 2*cx*cy*dy*ex + 2*cx*cy*dx*ey - pow(cx,2)*dy*ey + pow(cy,2)*dy*ey);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,aw*bw*dw*ew);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,cw*cw*aw*bw);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,cw*cw*dw*ew);

  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][j]->Add(qcnt_,aw*aw);
  qanal[anal].wnB[bin][j]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][j]->Fill(0.,ax*bx*dx*ex - ay*by*dx*ex + ay*bx*dy*ex + ax*by*dy*ex + ay*bx*dx*ey + ax*by*dx*ey - ax*bx*dy*ey + ay*by*dy*ey);
  qanal[anal].qCA[bin][j]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2));
  qanal[anal].qCB[bin][j]->Fill(0.,pow(cx,2)*dx*ex - pow(cy,2)*dx*ex + 2*cx*cy*dy*ex + 2*cx*cy*dx*ey - pow(cx,2)*dy*ey + pow(cy,2)*dy*ey);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,aw*bw*dw*ew);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,cw*cw*aw*bw);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,cw*cw*dw*ew);

}
void Fill_N42C(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp2];
  double ay = qy[HFp2];
  double aw = sumw[HFp2];

  double bx = qx[HFm2];
  double by = qy[HFm2];
  double bw = sumw[HFm2];

  double cx = qx[trackmid2];
  double cy = qy[trackmid2];
  double cw = sumw[trackmid2];

  double dx = qx[trackp218];
  double dy = qy[trackp218];
  double dw = sumw[trackp218];

  double ex = qx[trackm218];
  double ey = qy[trackm218];
  double ew = sumw[trackm218];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;
  if(pow(dx,2)+pow(dy,2) < 1e-6) return;
  if(pow(ex,2)+pow(ey,2) < 1e-6) return;


  qanal[anal].qA[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][0]->Add(qcnt_,aw*bw);
  qanal[anal].wnB[bin][0]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][0]->Fill(0.,ax*bx*dx*ex - ay*by*dx*ex + ay*bx*dy*ex + ax*by*dy*ex + ay*bx*dx*ey + ax*by*dx*ey - ax*bx*dy*ey + ay*by*dy*ey);
  qanal[anal].qCA[bin][0]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2));
  qanal[anal].qCB[bin][0]->Fill(0.,pow(cx,2)*dx*ex - pow(cy,2)*dx*ex + 2*cx*cy*dy*ex + 2*cx*cy*dx*ey - pow(cx,2)*dy*ey + pow(cy,2)*dy*ey);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,aw*bw*dw*ew);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,cw*cw*aw*bw);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,cw*cw*dw*ew);

  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][j]->Add(qcnt_,aw*aw);
  qanal[anal].wnB[bin][j]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][j]->Fill(0.,ax*bx*dx*ex - ay*by*dx*ex + ay*bx*dy*ex + ax*by*dy*ex + ay*bx*dx*ey + ax*by*dx*ey - ax*bx*dy*ey + ay*by*dy*ey);
  qanal[anal].qCA[bin][j]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2));
  qanal[anal].qCB[bin][j]->Fill(0.,pow(cx,2)*dx*ex - pow(cy,2)*dx*ex + 2*cx*cy*dy*ex + 2*cx*cy*dx*ey - pow(cx,2)*dy*ey + pow(cy,2)*dy*ey);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,aw*bw*dw*ew);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,cw*cw*aw*bw);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,cw*cw*dw*ew);

}

void Fill_N523(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFp3];
  double a3y = qy[HFp3];
  double a3w = sumw[HFp3];

  double b2x = qx[HFm2];
  double b2y = qy[HFm2];
  double b2w = sumw[HFm2];
  double b3x = qx[HFm3];
  double b3y = qy[HFm3];
  double b3w = sumw[HFm3];

  double c2x = qx[trackmid2];
  double c2y = qy[trackmid2];
  double c2w = sumw[trackmid2];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];

  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;


  qanal[anal].qA[bin][0]->Add(qxtrk_,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][0]->Add(qytrk_,a2y*a3x+a2x*a3y);
  qanal[anal].qB[bin][0]->Add(qxtrk_,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][0]->Add(qytrk_,b2y*b3x+b2x*b3y);
  qanal[anal].wnA[bin][0]->Add(qcnt_,a2w*a3w);
  qanal[anal].wnB[bin][0]->Add(qcnt_,b2w*b3w);
  qanal[anal].qBA[bin][0]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c3w*b2w*b3w);
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][j]->Add(qytrk_,a2y*a3x+a2x*a3y);
  qanal[anal].qB[bin][j]->Add(qxtrk_,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][j]->Add(qytrk_,b2y*b3x+b2x*b3y);
  qanal[anal].wnA[bin][j]->Add(qcnt_,a2w*a3w);
  qanal[anal].wnB[bin][j]->Add(qcnt_,b2w*b3w);
  qanal[anal].qBA[bin][j]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][j]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][j]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,c2w*c3w*b2w*b3w);
}
void Fill_N523A(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFm3];
  double a3y = qy[HFm3];
  double a3w = sumw[HFm3];

  double b2x = qx[HFm2];
  double b2y = qy[HFm2];
  double b2w = sumw[HFm2];
  double b3x = qx[HFp3];
  double b3y = qy[HFp3];
  double b3w = sumw[HFp3];

  double c2x = qx[trackmid2];
  double c2y = qy[trackmid2];
  double c2w = sumw[trackmid2];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];


  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][0]->Add(qytrk_,a2y*a3x+a2x*a3y);
  qanal[anal].qB[bin][0]->Add(qxtrk_,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][0]->Add(qytrk_,b2y*b3x+b2x*b3y);
  qanal[anal].wnA[bin][0]->Add(qcnt_,a2w*a3w);
  qanal[anal].wnB[bin][0]->Add(qcnt_,b2w*b3w);
  qanal[anal].qBA[bin][0]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c3w*b2w*b3w);
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][j]->Add(qytrk_,a2y*a3x+a2x*a3y);
  qanal[anal].qB[bin][j]->Add(qxtrk_,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][j]->Add(qytrk_,b2y*b3x+b2x*b3y);
  qanal[anal].wnA[bin][j]->Add(qcnt_,a2w*a3w);
  qanal[anal].wnB[bin][j]->Add(qcnt_,b2w*b3w);
  qanal[anal].qBA[bin][j]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][j]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][j]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,c2w*c3w*b2w*b3w);
}


void Fill_N63(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double Ax = qx[HFp3];
  double Ay = qy[HFp3];
  double wA = sumw[HFp3];

  double Bx = qx[HFm3];
  double By = qy[HFm3];
  double wB = sumw[HFm3];

  double Cx = qx[trackmid3];
  double Cy = qy[trackmid3];
  double wC = sumw[trackmid3];

  if(pow(Ax,2)+pow(Ay,2) < 1e-6) return;
  if(pow(Bx,2)+pow(By,2) < 1e-6) return;
  if(pow(Cx,2)+pow(Cy,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,pow(Ax,2)-pow(Ay,2));
  qanal[anal].qA[bin][0]->Add(qytrk_,2*Ax*Ay);
  qanal[anal].qB[bin][0]->Add(qxtrk_,pow(Bx,2)-pow(By,2));
  qanal[anal].qB[bin][0]->Add(qytrk_,2*Bx*By);
  qanal[anal].wnA[bin][0]->Add(qcnt_,pow(wA,2));
  qanal[anal].wnB[bin][0]->Add(qcnt_,pow(wB,2));
  qanal[anal].qBA[bin][0]->Fill(0.,pow(Ax,2)*(pow(Bx,2)-pow(By,2)) + pow(Ay,2)*(pow(By,2)-pow(Bx,2))+4*Ax*Ay*Bx*By);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(Ax,2)*(pow(Cx,2)-pow(Cy,2)) + pow(Ay,2)*(pow(Cy,2)-pow(Cx,2))+4*Ax*Ay*Cx*Cy);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(Bx,2)*(pow(Cx,2)-pow(Cy,2)) + pow(By,2)*(pow(Cy,2)-pow(Cx,2))+4*Bx*By*Cx*Cy);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wA*wB,2));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wA*wC,2));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wB*wC,2));
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,pow(Ax,2)-pow(Ay,2));
  qanal[anal].qA[bin][j]->Add(qytrk_,2*Ax*Ay);
  qanal[anal].qB[bin][j]->Add(qxtrk_,pow(Bx,2)-pow(By,2));
  qanal[anal].qB[bin][j]->Add(qytrk_,2*Bx*By);
  qanal[anal].wnA[bin][j]->Add(qcnt_,pow(wA,2));
  qanal[anal].wnB[bin][j]->Add(qcnt_,pow(wB,2));
  qanal[anal].qBA[bin][j]->Fill(0.,pow(Ax,2)*(pow(Bx,2)-pow(By,2)) + pow(Ay,2)*(pow(By,2)-pow(Bx,2))+4*Ax*Ay*Bx*By);
  qanal[anal].qCA[bin][j]->Fill(0.,pow(Ax,2)*(pow(Cx,2)-pow(Cy,2)) + pow(Ay,2)*(pow(Cy,2)-pow(Cx,2))+4*Ax*Ay*Cx*Cy);
  qanal[anal].qCB[bin][j]->Fill(0.,pow(Bx,2)*(pow(Cx,2)-pow(Cy,2)) + pow(By,2)*(pow(Cy,2)-pow(Cx,2))+4*Bx*By*Cx*Cy);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wA*wB,2));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wA*wC,2));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wB*wC,2));
}

void Fill_N63A(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp3];
  double ay = qy[HFp3];
  double aw = sumw[HFp3];

  double bx = qx[HFm3];
  double by = qy[HFm3];
  double bw = sumw[HFm3];

  double cx = qx[trackmid3];
  double cy = qy[trackmid3];
  double cw = sumw[trackmid3];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][0]->Add(qcnt_,aw*bw);
  qanal[anal].wnB[bin][0]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(aw*bw,2));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(aw*cw,2));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(bw*cw,2));

  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][j]->Add(qcnt_,aw*aw);
  qanal[anal].wnB[bin][j]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by);
  qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy);
  qanal[anal].qCB[bin][j]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(aw*bw,2));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(aw*cw,2));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(bw*cw,2));

}
void Fill_N63B(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp3];
  double ay = qy[HFp3];
  double aw = sumw[HFp3];

  double bx = qx[HFm3];
  double by = qy[HFm3];
  double bw = sumw[HFm3];

  double cx = qx[trackmid3];
  double cy = qy[trackmid3];
  double cw = sumw[trackmid3];

  double dx = qx[trackp310];
  double dy = qy[trackp310];
  double dw = sumw[trackp310];

  double ex = qx[trackm310];
  double ey = qy[trackm310];
  double ew = sumw[trackm310];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;
  if(pow(dx,2)+pow(dy,2) < 1e-6) return;
  if(pow(ex,2)+pow(ey,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][0]->Add(qcnt_,aw*bw);
  qanal[anal].wnB[bin][0]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][0]->Fill(0.,ax*bx*dx*ex - ay*by*dx*ex + ay*bx*dy*ex + ax*by*dy*ex + ay*bx*dx*ey + ax*by*dx*ey - ax*bx*dy*ey + ay*by*dy*ey);
  qanal[anal].qCA[bin][0]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2));
  qanal[anal].qCB[bin][0]->Fill(0.,pow(cx,2)*dx*ex - pow(cy,2)*dx*ex + 2*cx*cy*dy*ex + 2*cx*cy*dx*ey - pow(cx,2)*dy*ey + pow(cy,2)*dy*ey);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,aw*bw*dw*ew);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,cw*cw*aw*bw);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,cw*cw*dw*ew);

  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][j]->Add(qcnt_,aw*aw);
  qanal[anal].wnB[bin][j]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][j]->Fill(0.,ax*bx*dx*ex - ay*by*dx*ex + ay*bx*dy*ex + ax*by*dy*ex + ay*bx*dx*ey + ax*by*dx*ey - ax*bx*dy*ey + ay*by*dy*ey);
  qanal[anal].qCA[bin][j]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2));
  qanal[anal].qCB[bin][j]->Fill(0.,pow(cx,2)*dx*ex - pow(cy,2)*dx*ex + 2*cx*cy*dy*ex + 2*cx*cy*dx*ey - pow(cx,2)*dy*ey + pow(cy,2)*dy*ey);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,aw*bw*dw*ew);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,cw*cw*aw*bw);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,cw*cw*dw*ew);

}
void Fill_N63C(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp3];
  double ay = qy[HFp3];
  double aw = sumw[HFp3];

  double bx = qx[HFm3];
  double by = qy[HFm3];
  double bw = sumw[HFm3];

  double cx = qx[trackmid3];
  double cy = qy[trackmid3];
  double cw = sumw[trackmid3];

  double dx = qx[trackp318];
  double dy = qy[trackp318];
  double dw = sumw[trackp318];

  double ex = qx[trackm318];
  double ey = qy[trackm318];
  double ew = sumw[trackm318];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;
  if(pow(dx,2)+pow(dy,2) < 1e-6) return;
  if(pow(ex,2)+pow(ey,2) < 1e-6) return;


  qanal[anal].qA[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][0]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][0]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][0]->Add(qcnt_,aw*bw);
  qanal[anal].wnB[bin][0]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][0]->Fill(0.,ax*bx*dx*ex - ay*by*dx*ex + ay*bx*dy*ex + ax*by*dy*ex + ay*bx*dx*ey + ax*by*dx*ey - ax*bx*dy*ey + ay*by*dy*ey);
  qanal[anal].qCA[bin][0]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2));
  qanal[anal].qCB[bin][0]->Fill(0.,pow(cx,2)*dx*ex - pow(cy,2)*dx*ex + 2*cx*cy*dy*ex + 2*cx*cy*dx*ey - pow(cx,2)*dy*ey + pow(cy,2)*dy*ey);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,aw*bw*dw*ew);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,cw*cw*aw*bw);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,cw*cw*dw*ew);

  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qA[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].qB[bin][j]->Add(qxtrk_,ax*bx-ay*by);
  qanal[anal].qB[bin][j]->Add(qytrk_,ay*bx+ax*by);
  qanal[anal].wnA[bin][j]->Add(qcnt_,aw*aw);
  qanal[anal].wnB[bin][j]->Add(qcnt_,aw*bw);

  qanal[anal].qBA[bin][j]->Fill(0.,ax*bx*dx*ex - ay*by*dx*ex + ay*bx*dy*ex + ax*by*dy*ex + ay*bx*dx*ey + ax*by*dx*ey - ax*bx*dy*ey + ay*by*dy*ey);
  qanal[anal].qCA[bin][j]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2));
  qanal[anal].qCB[bin][j]->Fill(0.,pow(cx,2)*dx*ex - pow(cy,2)*dx*ex + 2*cx*cy*dy*ex + 2*cx*cy*dx*ey - pow(cx,2)*dy*ey + pow(cy,2)*dy*ey);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,aw*bw*dw*ew);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,cw*cw*aw*bw);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,cw*cw*dw*ew);

}

/* void Fill_N63(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){ */
/*   double ax = qx[HFp3]; */
/*   double ay = qy[HFp3]; */
/*   double wa = sumw[HFp3]; */

/*   double bx = qx[HFm3]; */
/*   double by = qy[HFm3]; */
/*   double wb = sumw[HFm3]; */

/*   double cx = qx[trackmid3]; */
/*   double cy = qy[trackmid3]; */
/*   double wc = sumw[trackmid3]; */


/*   if(pow(ax,2)+pow(ay,2) < 1e-6) return; */
/*   if(pow(bx,2)+pow(by,2) < 1e-6) return; */
/*   if(pow(cx,2)+pow(cy,2) < 1e-6) return; */

/*   qanal[anal].qA[bin][0]->Add(qxtrk_,pow(ax,2)-pow(ay,2)); */
/*   qanal[anal].qA[bin][0]->Add(qytrk_,2*ax*ay); */
/*   qanal[anal].qB[bin][0]->Add(qxtrk_,pow(bx,2)-pow(by,2)); */
/*   qanal[anal].qB[bin][0]->Add(qytrk_,2*bx*by); */
/*   qanal[anal].wnA[bin][0]->Add(qcnt_,pow(wa,2)); */
/*   qanal[anal].wnB[bin][0]->Add(qcnt_,pow(wb,2)); */
/*   qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by); */
/*   qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy); */
/*   qanal[anal].qCB[bin][0]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy); */
/*   qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wa*wb,2)); */
/*   qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wa*wc,2)); */
/*   qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wb*wc,2)); */
  
/*   int j=(int)(ran->Uniform(0,9.999))+1; */
/*   qanal[anal].qA[bin][j]->Add(qxtrk_,pow(ax,2)-pow(ay,2)); */
/*   qanal[anal].qA[bin][j]->Add(qytrk_,2*ax*ay); */
/*   qanal[anal].qB[bin][j]->Add(qxtrk_,pow(bx,2)-pow(by,2)); */
/*   qanal[anal].qB[bin][j]->Add(qytrk_,2*bx*by); */
/*   qanal[anal].wnA[bin][j]->Add(qcnt_,pow(wa,2)); */
/*   qanal[anal].wnB[bin][j]->Add(qcnt_,pow(wb,2)); */
/*   qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by); */
/*   qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy); */
/*   qanal[anal].qCB[bin][j]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy); */
/*   qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wa*wb,2)); */
/*   qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wa*wc,2)); */
/*   qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wb*wc,2)); */


/* } */




/* void Fill_N63A(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){ */
/*   double ax = qx[HFp3]; */
/*   double ay = qy[HFp3]; */
/*   double aw = sumw[HFp3]; */

/*   double bx = qx[HFm3]; */
/*   double by = qy[HFm3]; */
/*   double bw = sumw[HFm3]; */

/*   double cx = qx[trackmid3]; */
/*   double cy = qy[trackmid3]; */
/*   double cw = sumw[trackmid3]; */


/*   if(pow(ax,2)+pow(ay,2) < 1e-6) return; */
/*   if(pow(bx,2)+pow(by,2) < 1e-6) return; */
/*   if(pow(cx,2)+pow(cy,2) < 1e-6) return; */

/*   qanal[anal].qA[bin][0]->Add(qxtrk_,ax*bx-ay*by); */
/*   qanal[anal].qA[bin][0]->Add(qytrk_,ay*bx+ax*by); */
/*   qanal[anal].qB[bin][0]->Add(qxtrk_,ax*bx-ay*by); */
/*   qanal[anal].qB[bin][0]->Add(qytrk_,ay*bx+ax*by); */
/*   qanal[anal].wnA[bin][0]->Add(qcnt_,aw*bw); */
/*   qanal[anal].wnB[bin][0]->Add(qcnt_,aw*bw); */

/*   qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,2)*pow(bx,2) + pow(ay,2)*pow(bx,2) + pow(ax,2)*pow(by,2) + pow(ay,2)*pow(by,2)); */
/*   qanal[anal].qCA[bin][0]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2)); */
/*   qanal[anal].qCB[bin][0]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2)); */
/*   qanal[anal].qBAcnt[bin][0]->Fill(0.,bw*aw*bw*aw); */
/*   qanal[anal].qCAcnt[bin][0]->Fill(0.,cw*cw*aw*bw); */
/*   qanal[anal].qCBcnt[bin][0]->Fill(0.,cw*cw*bw*aw); */


  
/*   int j=(int)(ran->Uniform(0,9.999))+1; */
/*   qanal[anal].qA[bin][j]->Add(qxtrk_,ax*bx-ay*by); */
/*   qanal[anal].qA[bin][j]->Add(qytrk_,ay*bx+ax*by); */
/*   qanal[anal].qB[bin][j]->Add(qxtrk_,ax*bx-ay*by); */
/*   qanal[anal].qB[bin][j]->Add(qytrk_,ay*bx+ax*by); */
/*   qanal[anal].wnA[bin][j]->Add(qcnt_,aw*aw); */
/*   qanal[anal].wnB[bin][j]->Add(qcnt_,aw*bw); */

/*   qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,2)*pow(bx,2) + pow(ay,2)*pow(bx,2) + pow(ax,2)*pow(by,2) + pow(ay,2)*pow(by,2)); */
/*   qanal[anal].qCA[bin][j]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2)); */
/*   qanal[anal].qCB[bin][j]->Fill(0.,ax*bx*pow(cx,2) - ay*by*pow(cx,2) + 2*ay*bx*cx*cy + 2*ax*by*cx*cy - ax*bx*pow(cy,2) + ay*by*pow(cy,2)); */
/*   qanal[anal].qBAcnt[bin][j]->Fill(0.,bw*aw*bw*aw); */
/*   qanal[anal].qCAcnt[bin][j]->Fill(0.,cw*cw*aw*bw); */
/*   qanal[anal].qCBcnt[bin][j]->Fill(0.,cw*cw*bw*aw); */
/* } */






void Fill_N62(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp2];
  double ay = qy[HFp2];
  double wa = sumw[HFp2];

  double bx = qx[HFm2];
  double by = qy[HFm2];
  double wb = sumw[HFm2];

  double cx = qx[trackmid2];
  double cy = qy[trackmid2];
  double wc = sumw[trackmid2];


  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,pow(ax,3)-3*ax*pow(ay,2));
  qanal[anal].qA[bin][0]->Add(qytrk_,3*ay*pow(ax,2)-pow(ay,3));
  qanal[anal].qB[bin][0]->Add(qxtrk_,pow(bx,3)-3*bx*pow(by,2));
  qanal[anal].qB[bin][0]->Add(qytrk_,3*by*pow(bx,2)-pow(by,3));
  qanal[anal].wnA[bin][0]->Add(qcnt_,pow(wa,3));
  qanal[anal].wnB[bin][0]->Add(qcnt_,pow(wb,3));
  qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,3)*pow(bx,3) - 3*ax*pow(ay,2)*pow(bx,3) + 9*pow(ax,2)*ay*pow(bx,2)*by - 3*pow(ay,3)*pow(bx,2)*by - 3*pow(ax,3)*bx*pow(by,2) + 
   9*ax*pow(ay,2)*bx*pow(by,2) - 3*pow(ax,2)*ay*pow(by,3) + pow(ay,3)*pow(by,3));
  qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,3)*pow(cx,3) - 3*ax*pow(ay,2)*pow(cx,3) + 9*pow(ax,2)*ay*pow(cx,2)*cy - 3*pow(ay,3)*pow(cx,2)*cy - 3*pow(ax,3)*cx*pow(cy,2) + 
   9*ax*pow(ay,2)*cx*pow(cy,2) - 3*pow(ax,2)*ay*pow(cy,3) + pow(ay,3)*pow(cy,3));
  qanal[anal].qCB[bin][0]->Fill(0.,pow(bx,3)*pow(cx,3) - 3*bx*pow(by,2)*pow(cx,3) + 9*pow(bx,2)*by*pow(cx,2)*cy - 3*pow(by,3)*pow(cx,2)*cy - 3*pow(bx,3)*cx*pow(cy,2) + 
   9*bx*pow(by,2)*cx*pow(cy,2) - 3*pow(bx,2)*by*pow(cy,3) + pow(by,3)*pow(cy,3));
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wa*wb,3));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wa*wc,3));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wb*wc,3));
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,pow(ax,3)-3*ax*pow(ay,2));
  qanal[anal].qA[bin][j]->Add(qytrk_,3*ay*pow(ax,2)-pow(ay,3));
  qanal[anal].qB[bin][j]->Add(qxtrk_,pow(bx,3)-3*bx*pow(by,2));
  qanal[anal].qB[bin][j]->Add(qytrk_,3*by*pow(bx,2)-pow(by,3));
  qanal[anal].wnA[bin][j]->Add(qcnt_,pow(wa,3));
  qanal[anal].wnB[bin][j]->Add(qcnt_,pow(wb,3));
  qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,3)*pow(bx,3) - 3*ax*pow(ay,2)*pow(bx,3) + 9*pow(ax,2)*ay*pow(bx,2)*by - 3*pow(ay,3)*pow(bx,2)*by - 3*pow(ax,3)*bx*pow(by,2) + 
   9*ax*pow(ay,2)*bx*pow(by,2) - 3*pow(ax,2)*ay*pow(by,3) + pow(ay,3)*pow(by,3));
  qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,3)*pow(cx,3) - 3*ax*pow(ay,2)*pow(cx,3) + 9*pow(ax,2)*ay*pow(cx,2)*cy - 3*pow(ay,3)*pow(cx,2)*cy - 3*pow(ax,3)*cx*pow(cy,2) + 
   9*ax*pow(ay,2)*cx*pow(cy,2) - 3*pow(ax,2)*ay*pow(cy,3) + pow(ay,3)*pow(cy,3));
  qanal[anal].qCB[bin][j]->Fill(0.,pow(bx,3)*pow(cx,3) - 3*bx*pow(by,2)*pow(cx,3) + 9*pow(bx,2)*by*pow(cx,2)*cy - 3*pow(by,3)*pow(cx,2)*cy - 3*pow(bx,3)*cx*pow(cy,2) + 
   9*bx*pow(by,2)*cx*pow(cy,2) - 3*pow(bx,2)*by*pow(cy,3) + pow(by,3)*pow(cy,3));
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wa*wb,3));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wa*wc,3));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wb*wc,3));
 

}

void Fill_N62A(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp2];
  double ay = qy[HFp2];
  double wa = sumw[HFp2];

  double bx = qx[HFm2];
  double by = qy[HFm2];
  double wb = sumw[HFm2];

  double cx = qx[trackmid2];
  double cy = qy[trackmid2];
  double wc = sumw[trackmid2]; 


  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,ax*ax*bx-ay*ay*bx-2*ax*ay*by);
  qanal[anal].qA[bin][0]->Add(qytrk_,2*ax*ay*bx+ax*ax*by-ay*ay*by);
  qanal[anal].qB[bin][0]->Add(qxtrk_,bx*bx*ax-by*by*ax-2*bx*by*ay);
  qanal[anal].qB[bin][0]->Add(qytrk_,2*bx*by*ax+bx*bx*ay-by*by*ay);
  qanal[anal].wnA[bin][0]->Add(qcnt_,wa*wa*wb);
  qanal[anal].wnB[bin][0]->Add(qcnt_,wb*wb*wa);
  qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,3)*pow(bx,3) + ax*pow(ay,2)*pow(bx,3) + pow(ax,2)*ay*pow(bx,2)*by + pow(ay,3)*pow(bx,2)*by + pow(ax,3)*bx*pow(by,2) + 
   ax*pow(ay,2)*bx*pow(by,2) + pow(ax,2)*ay*pow(by,3) + pow(ay,3)*pow(by,3));
  qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,2)*bx*pow(cx,3) - pow(ay,2)*bx*pow(cx,3) - 2*ax*ay*by*pow(cx,3) + 6*ax*ay*bx*pow(cx,2)*cy + 3*pow(ax,2)*by*pow(cx,2)*cy - 
   3*pow(ay,2)*by*pow(cx,2)*cy - 3*pow(ax,2)*bx*cx*pow(cy,2) + 3*pow(ay,2)*bx*cx*pow(cy,2) + 6*ax*ay*by*cx*pow(cy,2) - 2*ax*ay*bx*pow(cy,3) - 
   pow(ax,2)*by*pow(cy,3) + pow(ay,2)*by*pow(cy,3));
  qanal[anal].qCB[bin][0]->Fill(0.,ax*pow(bx,2)*pow(cx,3) - 2*ay*bx*by*pow(cx,3) - ax*pow(by,2)*pow(cx,3) + 3*ay*pow(bx,2)*pow(cx,2)*cy + 6*ax*bx*by*pow(cx,2)*cy - 
   3*ay*pow(by,2)*pow(cx,2)*cy - 3*ax*pow(bx,2)*cx*pow(cy,2) + 6*ay*bx*by*cx*pow(cy,2) + 3*ax*pow(by,2)*cx*pow(cy,2) - 
   ay*pow(bx,2)*pow(cy,3) - 2*ax*bx*by*pow(cy,3) + ay*pow(by,2)*pow(cy,3));
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wb*wb*wa,2));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wc,3)*wa*wa*wb);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wc,3)*wb*wb*wa);
 
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,ax*ax*bx-ay*ay*bx-2*ax*ay*by);
  qanal[anal].qA[bin][j]->Add(qytrk_,2*ax*ay*bx+ax*ax*by-ay*ay*by);
  qanal[anal].qB[bin][j]->Add(qxtrk_,bx*bx*ax-by*by*ax-2*bx*by*ay);
  qanal[anal].qB[bin][j]->Add(qytrk_,2*bx*by*ax+bx*bx*ay-by*by*ay);
  qanal[anal].wnA[bin][j]->Add(qcnt_,wa*wa*wb);
  qanal[anal].wnB[bin][j]->Add(qcnt_,wb*wb*wa);
  qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,3)*pow(bx,3) + ax*pow(ay,2)*pow(bx,3) + pow(ax,2)*ay*pow(bx,2)*by + pow(ay,3)*pow(bx,2)*by + pow(ax,3)*bx*pow(by,2) + 
   ax*pow(ay,2)*bx*pow(by,2) + pow(ax,2)*ay*pow(by,3) + pow(ay,3)*pow(by,3));
  qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,2)*bx*pow(cx,3) - pow(ay,2)*bx*pow(cx,3) - 2*ax*ay*by*pow(cx,3) + 6*ax*ay*bx*pow(cx,2)*cy + 3*pow(ax,2)*by*pow(cx,2)*cy - 
   3*pow(ay,2)*by*pow(cx,2)*cy - 3*pow(ax,2)*bx*cx*pow(cy,2) + 3*pow(ay,2)*bx*cx*pow(cy,2) + 6*ax*ay*by*cx*pow(cy,2) - 2*ax*ay*bx*pow(cy,3) - 
   pow(ax,2)*by*pow(cy,3) + pow(ay,2)*by*pow(cy,3));
  qanal[anal].qCB[bin][j]->Fill(0.,ax*pow(bx,2)*pow(cx,3) - 2*ay*bx*by*pow(cx,3) - ax*pow(by,2)*pow(cx,3) + 3*ay*pow(bx,2)*pow(cx,2)*cy + 6*ax*bx*by*pow(cx,2)*cy - 
   3*ay*pow(by,2)*pow(cx,2)*cy - 3*ax*pow(bx,2)*cx*pow(cy,2) + 6*ay*bx*by*cx*pow(cy,2) + 3*ax*pow(by,2)*cx*pow(cy,2) - 
   ay*pow(bx,2)*pow(cy,3) - 2*ax*bx*by*pow(cy,3) + ay*pow(by,2)*pow(cy,3));
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wb*wb*wa,2));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wc,3)*wa*wa*wb);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wc,3)*wb*wb*wa);
 

}
void Fill_N723(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFp3];
  double a3y = qy[HFp3];
  double a3w = sumw[HFp3];

  double b2x = qx[HFm2];
  double b2y = qy[HFm2];
  double b2w = sumw[HFm2];
  double b3x = qx[HFm3];
  double b3y = qy[HFm3];
  double b3w = sumw[HFm3];

  double c2x = qx[trackmid2];
  double c2y = qy[trackmid2];
  double c2w = sumw[trackmid2];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];


  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,pow(a2x,2)*a3x - pow(a2y,2)*a3x - 2*a2x*a2y*a3y);
  qanal[anal].qA[bin][0]->Add(qytrk_,2*a2x*a2y*a3x  + pow(a2x,2)*a3y - pow(a2y,2)*a3y );
  qanal[anal].qB[bin][0]->Add(qxtrk_,pow(b2x,2)*b3x - pow(b2y,2)*b3x - 2*b2x*b2y*b3y);
  qanal[anal].qB[bin][0]->Add(qytrk_, 2*b2x*b2y*b3x + pow(b2x,2)*b3y - pow(b2y,2)*b3y );
  qanal[anal].wnA[bin][0]->Add(qcnt_,a2w*a2w*a3w);
  qanal[anal].wnB[bin][0]->Add(qcnt_,b2w*b2w*b3w);
  qanal[anal].qBA[bin][0]->Fill(0.,pow(a2x,2)*a3x*pow(b2x,2)*b3x - pow(a2y,2)*a3x*pow(b2x,2)*b3x - 2*a2x*a2y*a3y*pow(b2x,2)*b3x + 4*a2x*a2y*a3x*b2x*b2y*b3x + 2*pow(a2x,2)*a3y*b2x*b2y*b3x - 
		    2*pow(a2y,2)*a3y*b2x*b2y*b3x - pow(a2x,2)*a3x*pow(b2y,2)*b3x + pow(a2y,2)*a3x*pow(b2y,2)*b3x + 2*a2x*a2y*a3y*pow(b2y,2)*b3x + 2*a2x*a2y*a3x*pow(b2x,2)*b3y + 
		    pow(a2x,2)*a3y*pow(b2x,2)*b3y - pow(a2y,2)*a3y*pow(b2x,2)*b3y - 2*pow(a2x,2)*a3x*b2x*b2y*b3y + 2*pow(a2y,2)*a3x*b2x*b2y*b3y + 4*a2x*a2y*a3y*b2x*b2y*b3y - 
		    2*a2x*a2y*a3x*pow(b2y,2)*b3y - pow(a2x,2)*a3y*pow(b2y,2)*b3y + pow(a2y,2)*a3y*pow(b2y,2)*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(a2x,2)*a3x*pow(c2x,2)*c3x - pow(a2y,2)*a3x*pow(c2x,2)*c3x - 2*a2x*a2y*a3y*pow(c2x,2)*c3x + 4*a2x*a2y*a3x*c2x*c2y*c3x + 2*pow(a2x,2)*a3y*c2x*c2y*c3x - 
		    2*pow(a2y,2)*a3y*c2x*c2y*c3x - pow(a2x,2)*a3x*pow(c2y,2)*c3x + pow(a2y,2)*a3x*pow(c2y,2)*c3x + 2*a2x*a2y*a3y*pow(c2y,2)*c3x + 2*a2x*a2y*a3x*pow(c2x,2)*c3y + 
		    pow(a2x,2)*a3y*pow(c2x,2)*c3y - pow(a2y,2)*a3y*pow(c2x,2)*c3y - 2*pow(a2x,2)*a3x*c2x*c2y*c3y + 2*pow(a2y,2)*a3x*c2x*c2y*c3y + 4*a2x*a2y*a3y*c2x*c2y*c3y - 
		    2*a2x*a2y*a3x*pow(c2y,2)*c3y - pow(a2x,2)*a3y*pow(c2y,2)*c3y + pow(a2y,2)*a3y*pow(c2y,2)*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(b2x,2)*b3x*pow(c2x,2)*c3x - pow(b2y,2)*b3x*pow(c2x,2)*c3x - 2*b2x*b2y*b3y*pow(c2x,2)*c3x + 4*b2x*b2y*b3x*c2x*c2y*c3x + 2*pow(b2x,2)*b3y*c2x*c2y*c3x - 
		    2*pow(b2y,2)*b3y*c2x*c2y*c3x - pow(b2x,2)*b3x*pow(c2y,2)*c3x + pow(b2y,2)*b3x*pow(c2y,2)*c3x + 2*b2x*b2y*b3y*pow(c2y,2)*c3x + 2*b2x*b2y*b3x*pow(c2x,2)*c3y + 
		    pow(b2x,2)*b3y*pow(c2x,2)*c3y - pow(b2y,2)*b3y*pow(c2x,2)*c3y - 2*pow(b2x,2)*b3x*c2x*c2y*c3y + 2*pow(b2y,2)*b3x*c2x*c2y*c3y + 4*b2x*b2y*b3y*c2x*c2y*c3y - 
		    2*b2x*b2y*b3x*pow(c2y,2)*c3y - pow(b2x,2)*b3y*pow(c2y,2)*c3y + pow(b2y,2)*b3y*pow(c2y,2)*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b2w*b3w*a2w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c2w*c3w*a2w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c2w*c3w*b2w*b2w*b3w);
 
 int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,pow(a2x,2)*a3x - pow(a2y,2)*a3x- 2*a2x*a2y*a3y);
  qanal[anal].qA[bin][j]->Add(qytrk_,2*a2x*a2y*a3x + pow(a2x,2)*a3y- pow(a2y,2)*a3y );
  qanal[anal].qB[bin][j]->Add(qxtrk_,pow(b2x,2)*b3x - pow(b2y,2)*b3x - 2*b2x*b2y*b3y);
  qanal[anal].qB[bin][j]->Add(qytrk_, 2*b2x*b2y*b3x + pow(b2x,2)*b3y - pow(b2y,2)*b3y );
  qanal[anal].wnA[bin][j]->Add(qcnt_,a2w*a2w*a3w);
  qanal[anal].wnB[bin][j]->Add(qcnt_,b2w*a2w*b3w);
  qanal[anal].qBA[bin][j]->Fill(0.,pow(a2x,2)*a3x*pow(b2x,2)*b3x - pow(a2y,2)*a3x*pow(b2x,2)*b3x - 2*a2x*a2y*a3y*pow(b2x,2)*b3x + 4*a2x*a2y*a3x*b2x*b2y*b3x + 2*pow(a2x,2)*a3y*b2x*b2y*b3x - 
   2*pow(a2y,2)*a3y*b2x*b2y*b3x - pow(a2x,2)*a3x*pow(b2y,2)*b3x + pow(a2y,2)*a3x*pow(b2y,2)*b3x + 2*a2x*a2y*a3y*pow(b2y,2)*b3x + 2*a2x*a2y*a3x*pow(b2x,2)*b3y + 
   pow(a2x,2)*a3y*pow(b2x,2)*b3y - pow(a2y,2)*a3y*pow(b2x,2)*b3y - 2*pow(a2x,2)*a3x*b2x*b2y*b3y + 2*pow(a2y,2)*a3x*b2x*b2y*b3y + 4*a2x*a2y*a3y*b2x*b2y*b3y - 
   2*a2x*a2y*a3x*pow(b2y,2)*b3y - pow(a2x,2)*a3y*pow(b2y,2)*b3y + pow(a2y,2)*a3y*pow(b2y,2)*b3y);
  qanal[anal].qCA[bin][j]->Fill(0.,pow(a2x,2)*a3x*pow(c2x,2)*c3x - pow(a2y,2)*a3x*pow(c2x,2)*c3x - 2*a2x*a2y*a3y*pow(c2x,2)*c3x + 4*a2x*a2y*a3x*c2x*c2y*c3x + 2*pow(a2x,2)*a3y*c2x*c2y*c3x - 
   2*pow(a2y,2)*a3y*c2x*c2y*c3x - pow(a2x,2)*a3x*pow(c2y,2)*c3x + pow(a2y,2)*a3x*pow(c2y,2)*c3x + 2*a2x*a2y*a3y*pow(c2y,2)*c3x + 2*a2x*a2y*a3x*pow(c2x,2)*c3y + 
   pow(a2x,2)*a3y*pow(c2x,2)*c3y - pow(a2y,2)*a3y*pow(c2x,2)*c3y - 2*pow(a2x,2)*a3x*c2x*c2y*c3y + 2*pow(a2y,2)*a3x*c2x*c2y*c3y + 4*a2x*a2y*a3y*c2x*c2y*c3y - 
   2*a2x*a2y*a3x*pow(c2y,2)*c3y - pow(a2x,2)*a3y*pow(c2y,2)*c3y + pow(a2y,2)*a3y*pow(c2y,2)*c3y);
  qanal[anal].qCB[bin][j]->Fill(0.,pow(b2x,2)*b3x*pow(c2x,2)*c3x - pow(b2y,2)*b3x*pow(c2x,2)*c3x - 2*b2x*b2y*b3y*pow(c2x,2)*c3x + 4*b2x*b2y*b3x*c2x*c2y*c3x + 2*pow(b2x,2)*b3y*c2x*c2y*c3x - 
   2*pow(b2y,2)*b3y*c2x*c2y*c3x - pow(b2x,2)*b3x*pow(c2y,2)*c3x + pow(b2y,2)*b3x*pow(c2y,2)*c3x + 2*b2x*b2y*b3y*pow(c2y,2)*c3x + 2*b2x*b2y*b3x*pow(c2x,2)*c3y + 
   pow(b2x,2)*b3y*pow(c2x,2)*c3y - pow(b2y,2)*b3y*pow(c2x,2)*c3y - 2*pow(b2x,2)*b3x*c2x*c2y*c3y + 2*pow(b2y,2)*b3x*c2x*c2y*c3y + 4*b2x*b2y*b3y*c2x*c2y*c3y - 
   2*b2x*b2y*b3x*pow(c2y,2)*c3y - pow(b2x,2)*b3y*pow(c2y,2)*c3y + pow(b2y,2)*b3y*pow(c2y,2)*c3y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b2w*b2w*b3w*a2w*a2w*a3w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,c2w*c2w*c3w*a2w*a2w*a3w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,c2w*c2w*c3w*b2w*b2w*b3w);

}

void Fill_N723A(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFm3];
  double a3y = qy[HFm3];
  double a3w = sumw[HFm3];

  double b2x = qx[HFm2];
  double b2y = qy[HFm2];
  double b2w = sumw[HFm2];
  double b3x = qx[HFp3];
  double b3y = qy[HFp3];
  double b3w = sumw[HFp3];

  double c2x = qx[trackmid2];
  double c2y = qy[trackmid2];
  double c2w = sumw[trackmid2];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];

  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  qanal[anal].qA[bin][0]->Add(qxtrk_,pow(a2x,2)*a3x - pow(a2y,2)*a3x- 2*a2x*a2y*a3y);
  qanal[anal].qA[bin][0]->Add(qytrk_,2*a2x*a2y*a3x + pow(a2x,2)*a3y- pow(a2y,2)*a3y );
  qanal[anal].qB[bin][0]->Add(qxtrk_,pow(b2x,2)*b3x - pow(b2y,2)*b3x - 2*b2x*b2y*b3y);
  qanal[anal].qB[bin][0]->Add(qytrk_, 2*b2x*b2y*b3x + pow(b2x,2)*b3y - pow(b2y,2)*b3y );
  qanal[anal].wnA[bin][0]->Add(qcnt_,a2w*a2w*a3w);
  qanal[anal].wnB[bin][0]->Add(qcnt_,b2w*b2w*b3w);
  qanal[anal].qBA[bin][0]->Fill(0.,pow(a2x,2)*a3x*pow(b2x,2)*b3x - pow(a2y,2)*a3x*pow(b2x,2)*b3x - 2*a2x*a2y*a3y*pow(b2x,2)*b3x + 4*a2x*a2y*a3x*b2x*b2y*b3x + 2*pow(a2x,2)*a3y*b2x*b2y*b3x - 
   2*pow(a2y,2)*a3y*b2x*b2y*b3x - pow(a2x,2)*a3x*pow(b2y,2)*b3x + pow(a2y,2)*a3x*pow(b2y,2)*b3x + 2*a2x*a2y*a3y*pow(b2y,2)*b3x + 2*a2x*a2y*a3x*pow(b2x,2)*b3y + 
   pow(a2x,2)*a3y*pow(b2x,2)*b3y - pow(a2y,2)*a3y*pow(b2x,2)*b3y - 2*pow(a2x,2)*a3x*b2x*b2y*b3y + 2*pow(a2y,2)*a3x*b2x*b2y*b3y + 4*a2x*a2y*a3y*b2x*b2y*b3y - 
   2*a2x*a2y*a3x*pow(b2y,2)*b3y - pow(a2x,2)*a3y*pow(b2y,2)*b3y + pow(a2y,2)*a3y*pow(b2y,2)*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(a2x,2)*a3x*pow(c2x,2)*c3x - pow(a2y,2)*a3x*pow(c2x,2)*c3x - 2*a2x*a2y*a3y*pow(c2x,2)*c3x + 4*a2x*a2y*a3x*c2x*c2y*c3x + 2*pow(a2x,2)*a3y*c2x*c2y*c3x - 
   2*pow(a2y,2)*a3y*c2x*c2y*c3x - pow(a2x,2)*a3x*pow(c2y,2)*c3x + pow(a2y,2)*a3x*pow(c2y,2)*c3x + 2*a2x*a2y*a3y*pow(c2y,2)*c3x + 2*a2x*a2y*a3x*pow(c2x,2)*c3y + 
   pow(a2x,2)*a3y*pow(c2x,2)*c3y - pow(a2y,2)*a3y*pow(c2x,2)*c3y - 2*pow(a2x,2)*a3x*c2x*c2y*c3y + 2*pow(a2y,2)*a3x*c2x*c2y*c3y + 4*a2x*a2y*a3y*c2x*c2y*c3y - 
   2*a2x*a2y*a3x*pow(c2y,2)*c3y - pow(a2x,2)*a3y*pow(c2y,2)*c3y + pow(a2y,2)*a3y*pow(c2y,2)*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(b2x,2)*b3x*pow(c2x,2)*c3x - pow(b2y,2)*b3x*pow(c2x,2)*c3x - 2*b2x*b2y*b3y*pow(c2x,2)*c3x + 4*b2x*b2y*b3x*c2x*c2y*c3x + 2*pow(b2x,2)*b3y*c2x*c2y*c3x - 
   2*pow(b2y,2)*b3y*c2x*c2y*c3x - pow(b2x,2)*b3x*pow(c2y,2)*c3x + pow(b2y,2)*b3x*pow(c2y,2)*c3x + 2*b2x*b2y*b3y*pow(c2y,2)*c3x + 2*b2x*b2y*b3x*pow(c2x,2)*c3y + 
   pow(b2x,2)*b3y*pow(c2x,2)*c3y - pow(b2y,2)*b3y*pow(c2x,2)*c3y - 2*pow(b2x,2)*b3x*c2x*c2y*c3y + 2*pow(b2y,2)*b3x*c2x*c2y*c3y + 4*b2x*b2y*b3y*c2x*c2y*c3y - 
   2*b2x*b2y*b3x*pow(c2y,2)*c3y - pow(b2x,2)*b3y*pow(c2y,2)*c3y + pow(b2y,2)*b3y*pow(c2y,2)*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b2w*b3w*a2w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c2w*c3w*a2w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c2w*c3w*b2w*b2w*b3w);
 
 int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,pow(a2x,2)*a3x - pow(a2y,2)*a3x- 2*a2x*a2y*a3y);
  qanal[anal].qA[bin][j]->Add(qytrk_,2*a2x*a2y*a3x + pow(a2x,2)*a3y- pow(a2y,2)*a3y );
  qanal[anal].qB[bin][j]->Add(qxtrk_,pow(b2x,2)*b3x - pow(b2y,2)*b3x - 2*b2x*b2y*b3y);
  qanal[anal].qB[bin][j]->Add(qytrk_, 2*b2x*b2y*b3x + pow(b2x,2)*b3y - pow(b2y,2)*b3y );
  qanal[anal].wnA[bin][j]->Add(qcnt_,a2w*a2w*a3w);
  qanal[anal].wnB[bin][j]->Add(qcnt_,b2w*a2w*b3w);
  qanal[anal].qBA[bin][j]->Fill(0.,pow(a2x,2)*a3x*pow(b2x,2)*b3x - pow(a2y,2)*a3x*pow(b2x,2)*b3x - 2*a2x*a2y*a3y*pow(b2x,2)*b3x + 4*a2x*a2y*a3x*b2x*b2y*b3x + 2*pow(a2x,2)*a3y*b2x*b2y*b3x - 
   2*pow(a2y,2)*a3y*b2x*b2y*b3x - pow(a2x,2)*a3x*pow(b2y,2)*b3x + pow(a2y,2)*a3x*pow(b2y,2)*b3x + 2*a2x*a2y*a3y*pow(b2y,2)*b3x + 2*a2x*a2y*a3x*pow(b2x,2)*b3y + 
   pow(a2x,2)*a3y*pow(b2x,2)*b3y - pow(a2y,2)*a3y*pow(b2x,2)*b3y - 2*pow(a2x,2)*a3x*b2x*b2y*b3y + 2*pow(a2y,2)*a3x*b2x*b2y*b3y + 4*a2x*a2y*a3y*b2x*b2y*b3y - 
   2*a2x*a2y*a3x*pow(b2y,2)*b3y - pow(a2x,2)*a3y*pow(b2y,2)*b3y + pow(a2y,2)*a3y*pow(b2y,2)*b3y);
  qanal[anal].qCA[bin][j]->Fill(0.,pow(a2x,2)*a3x*pow(c2x,2)*c3x - pow(a2y,2)*a3x*pow(c2x,2)*c3x - 2*a2x*a2y*a3y*pow(c2x,2)*c3x + 4*a2x*a2y*a3x*c2x*c2y*c3x + 2*pow(a2x,2)*a3y*c2x*c2y*c3x - 
   2*pow(a2y,2)*a3y*c2x*c2y*c3x - pow(a2x,2)*a3x*pow(c2y,2)*c3x + pow(a2y,2)*a3x*pow(c2y,2)*c3x + 2*a2x*a2y*a3y*pow(c2y,2)*c3x + 2*a2x*a2y*a3x*pow(c2x,2)*c3y + 
   pow(a2x,2)*a3y*pow(c2x,2)*c3y - pow(a2y,2)*a3y*pow(c2x,2)*c3y - 2*pow(a2x,2)*a3x*c2x*c2y*c3y + 2*pow(a2y,2)*a3x*c2x*c2y*c3y + 4*a2x*a2y*a3y*c2x*c2y*c3y - 
   2*a2x*a2y*a3x*pow(c2y,2)*c3y - pow(a2x,2)*a3y*pow(c2y,2)*c3y + pow(a2y,2)*a3y*pow(c2y,2)*c3y);
  qanal[anal].qCB[bin][j]->Fill(0.,pow(b2x,2)*b3x*pow(c2x,2)*c3x - pow(b2y,2)*b3x*pow(c2x,2)*c3x - 2*b2x*b2y*b3y*pow(c2x,2)*c3x + 4*b2x*b2y*b3x*c2x*c2y*c3x + 2*pow(b2x,2)*b3y*c2x*c2y*c3x - 
   2*pow(b2y,2)*b3y*c2x*c2y*c3x - pow(b2x,2)*b3x*pow(c2y,2)*c3x + pow(b2y,2)*b3x*pow(c2y,2)*c3x + 2*b2x*b2y*b3y*pow(c2y,2)*c3x + 2*b2x*b2y*b3x*pow(c2x,2)*c3y + 
   pow(b2x,2)*b3y*pow(c2x,2)*c3y - pow(b2y,2)*b3y*pow(c2x,2)*c3y - 2*pow(b2x,2)*b3x*c2x*c2y*c3y + 2*pow(b2y,2)*b3x*c2x*c2y*c3y + 4*b2x*b2y*b3y*c2x*c2y*c3y - 
   2*b2x*b2y*b3x*pow(c2y,2)*c3y - pow(b2x,2)*b3y*pow(c2y,2)*c3y + pow(b2y,2)*b3y*pow(c2y,2)*c3y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b2w*b2w*b3w*a2w*a2w*a3w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,c2w*c2w*c3w*a2w*a2w*a3w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,c2w*c2w*c3w*b2w*b2w*b3w);

}
void Fill_D24(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp2];
  double ay = qy[HFp2];
  double wa = sumw[HFp2];

  double bx = qx[HFm2];
  double by = qy[HFm2];
  double wb = sumw[HFm2];

  double cx = qx[trackmid2];
  double cy = qy[trackmid2];
  double wc = sumw[trackmid2];


  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;

  qxy->Reset();
  qxy->Add(qxtrk_,1.);
  qxy->Multiply(qytrk_);

  qxt->Reset();
  qxt->Add(qxtrk_,1.);
  qxt->Multiply(qxtrk_);
  qyt->Reset();
  qyt->Add(qytrk_,1.);
  qyt->Multiply(qytrk_);
  qct->Reset();
  qct->Add(qcnt_,1.);
  qct->Multiply(qcnt_);

  qanal[anal].qA[bin][0]->Add(qxt, pow(ax,2) - pow(ay,2));
  qanal[anal].qA[bin][0]->Add(qyt, pow(ay,2) - pow(ax,2));
  qanal[anal].qA[bin][0]->Add(qxy, 4*ax*ay);
  
  qanal[anal].qB[bin][0]->Add(qxt, pow(bx,2) - pow(by,2));
  qanal[anal].qB[bin][0]->Add(qyt, pow(by,2) - pow(bx,2));
  qanal[anal].qB[bin][0]->Add(qxy, 4*bx*by);
  
  qanal[anal].wnA[bin][0]->Add(qct,pow(wa,2));
  qanal[anal].wnB[bin][0]->Add(qct,pow(wb,2));
  qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wa*wb,2));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wa*wc,2));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wb*wc,2));
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxt, pow(ax,2) - pow(ay,2));
  qanal[anal].qA[bin][j]->Add(qyt, pow(ay,2) - pow(ax,2));
  qanal[anal].qA[bin][j]->Add(qxy, 4*ax*ay);
  
  qanal[anal].qB[bin][j]->Add(qxt, pow(bx,2) - pow(by,2));
  qanal[anal].qB[bin][j]->Add(qyt, pow(by,2) - pow(bx,2));
  qanal[anal].qB[bin][j]->Add(qxy, 4*bx*by);
  
  qanal[anal].wnA[bin][j]->Add(qct,pow(wa,2));
  qanal[anal].wnB[bin][j]->Add(qct,pow(wb,2));
  qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by);
  qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy);
  qanal[anal].qCB[bin][j]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wa*wb,2));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wa*wc,2));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wb*wc,2));
  
}

void Fill_D24A(int anal, int bin, TH2D * qxtrk2_, TH2D * qytrk2_, TH2D * qxtrk3_, TH2D * qytrk3_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFm2];
  double a3y = qy[HFm2];
  double a3w = sumw[HFm2];
  
  double b2x = qx[HFm2];
  double b2y = qy[HFm2];
  double b2w = sumw[HFm2];
  double b3x = qx[HFp2];
  double b3y = qy[HFp2];
  double b3w = sumw[HFp2];
  
  double c2x = qx[trackmid2];
  double c2y = qy[trackmid2];
  double c2w = sumw[trackmid2];
  double c3x = qx[trackmid2];
  double c3y = qy[trackmid2];
  double c3w = sumw[trackmid2];


  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  qx2y3->Reset();
  qx2y3->Add(qxtrk2_,1.);
  qx2y3->Multiply(qytrk3_);

  qy2x3->Reset();
  qy2x3->Add(qytrk2_,1.);
  qy2x3->Multiply(qxtrk3_);


  qxt->Reset();
  qxt->Add(qxtrk2_,1.);
  qxt->Multiply(qxtrk3_);
  qyt->Reset();
  qyt->Add(qytrk2_,1.);
  qyt->Multiply(qytrk3_);
  qct->Reset();
  qct->Add(qcnt_,1.);
  qct->Multiply(qcnt_);

  qanal[anal].qA[bin][0]->Add(qxt,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][0]->Add(qyt,a2y*a3y-a2x*a3x);
  qanal[anal].qA[bin][0]->Add(qy2x3,  a2y*a3x+a2x*a3y);
  qanal[anal].qA[bin][0]->Add(qx2y3,  a2y*a3x+a2x*a3y);

  qanal[anal].qB[bin][0]->Add(qxt,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][0]->Add(qyt,b2y*b3y-b2x*b3x);
  qanal[anal].qB[bin][0]->Add(qy2x3,  b2y*b3x+b2x*b3y);
  qanal[anal].qB[bin][0]->Add(qx2y3,  b2y*b3x+b2x*b3y);

  qanal[anal].wnA[bin][0]->Add(qct,a2w*a3w);
  qanal[anal].wnB[bin][0]->Add(qct,b2w*b3w);

  qanal[anal].qBA[bin][0]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c3w*b2w*b3w);
  
  int j=(int)(ran->Uniform(0,9.999))+1;

  qanal[anal].qA[bin][j]->Add(qxt,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][j]->Add(qyt,a2y*a3y-a2x*a3x);
  qanal[anal].qA[bin][j]->Add(qy2x3,  a2y*a3x+a2x*a3y);
  qanal[anal].qA[bin][j]->Add(qx2y3,  a2y*a3x+a2x*a3y);

  qanal[anal].qB[bin][j]->Add(qxt,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][j]->Add(qyt,b2y*b3y-b2x*b3x);
  qanal[anal].qB[bin][j]->Add(qy2x3,  b2y*b3x+b2x*b3y);
  qanal[anal].qB[bin][j]->Add(qx2y3,  b2y*b3x+b2x*b3y);

  qanal[anal].wnA[bin][j]->Add(qct,a2w*a3w);
  qanal[anal].wnB[bin][j]->Add(qct,b2w*b3w);

  qanal[anal].qBA[bin][j]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][j]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][j]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,c2w*c3w*b2w*b3w);
  
  
}


void Fill_D26(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp2];
  double ay = qy[HFp2];
  double wa = sumw[HFp2];

  double bx = qx[HFm2];
  double by = qy[HFm2];
  double wb = sumw[HFm2];

  double cx = qx[trackmid2];
  double cy = qy[trackmid2];
  double wc = sumw[trackmid2];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;

  qxxy->Reset();
  qxxy->Add(qxtrk_,1.);
  qxxy->Multiply(qxtrk_);
  qxxy->Multiply(qytrk_);

  qxyy->Reset();
  qxyy->Add(qxtrk_,1.);
  qxyy->Multiply(qytrk_);
  qxyy->Multiply(qytrk_);

  qxt->Reset();
  qxt->Add(qxtrk_,1.);
  qxt->Multiply(qxtrk_);
  qxt->Multiply(qxtrk_);
  qyt->Reset();
  qyt->Add(qytrk_,1.);
  qyt->Multiply(qytrk_);
  qyt->Multiply(qytrk_);
  qct->Reset();
  qct->Add(qcnt_,1.);
  qct->Multiply(qcnt_);
  qct->Multiply(qcnt_);


  qanal[anal].qA[bin][0]->Add(qxt, pow(ax,3) - 3*ax*pow(ay,2));
  qanal[anal].qA[bin][0]->Add(qyt, pow(ay,3) - 3*ay*pow(ax,2));
  qanal[anal].qA[bin][0]->Add(qxxy, 9*pow(ax,2)*ay - 3*pow(ay,3));
  qanal[anal].qA[bin][0]->Add(qxyy, 9*pow(ay,2)*ax - 3*pow(ax,3));

  qanal[anal].qB[bin][0]->Add(qxt, pow(bx,3) - 3*bx*pow(by,2));
  qanal[anal].qB[bin][0]->Add(qyt, pow(by,3) - 3*by*pow(bx,2));
  qanal[anal].qB[bin][0]->Add(qxxy, 9*pow(bx,2)*by - 3*pow(by,3));
  qanal[anal].qB[bin][0]->Add(qxyy, 9*pow(by,2)*bx - 3*pow(bx,3));
  
  qanal[anal].wnA[bin][0]->Add(qct,pow(wa,3));
  qanal[anal].wnB[bin][0]->Add(qct,pow(wb,3));
  qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,3)*pow(bx,3) - 2*pow(ax,2)*ay*pow(bx,3) - ax*pow(ay,2)*pow(bx,3) + 3*pow(ax,3)*pow(bx,2)*by + 6*pow(ax,2)*ay*pow(bx,2)*by - 
		    3*ax*pow(ay,2)*pow(bx,2)*by - 3*pow(ax,3)*bx*pow(by,2) + 6*pow(ax,2)*ay*bx*pow(by,2) + 3*ax*pow(ay,2)*bx*pow(by,2) - 
		    pow(ax,3)*pow(by,3) - 2*pow(ax,2)*ay*pow(by,3) + ax*pow(ay,2)*pow(by,3));
  qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,3)*pow(cx,3) - 2*pow(ax,2)*ay*pow(cx,3) - ax*pow(ay,2)*pow(cx,3) + 3*pow(ax,3)*pow(cx,2)*cy + 6*pow(ax,2)*ay*pow(cx,2)*cy - 
		    3*ax*pow(ay,2)*pow(cx,2)*cy - 3*pow(ax,3)*cx*pow(cy,2) + 6*pow(ax,2)*ay*cx*pow(cy,2) + 3*ax*pow(ay,2)*cx*pow(cy,2) - 
		    pow(ax,3)*pow(cy,3) - 2*pow(ax,2)*ay*pow(cy,3) + ax*pow(ay,2)*pow(cy,3));
  qanal[anal].qCB[bin][0]->Fill(0.,pow(bx,3)*pow(cx,3) - 2*pow(bx,2)*by*pow(cx,3) - bx*pow(by,2)*pow(cx,3) + 3*pow(bx,3)*pow(cx,2)*cy + 6*pow(bx,2)*by*pow(cx,2)*cy - 
		    3*bx*pow(by,2)*pow(cx,2)*cy - 3*pow(bx,3)*cx*pow(cy,2) + 6*pow(bx,2)*by*cx*pow(cy,2) + 3*bx*pow(by,2)*cx*pow(cy,2) - 
		    pow(bx,3)*pow(cy,3) - 2*pow(bx,2)*by*pow(cy,3) + bx*pow(by,2)*pow(cy,3));
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wa*wb,3));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wa*wc,3));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wb*wc,3));
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxt, pow(ax,3) - 3*ax*pow(ay,2));
  qanal[anal].qA[bin][j]->Add(qyt, pow(ay,3) - 3*ay*pow(ax,2));
  qanal[anal].qA[bin][j]->Add(qxxy, 9*pow(ax,2)*ay - 3*pow(ay,3));
  qanal[anal].qA[bin][j]->Add(qxyy, 9*pow(ay,2)*ax - 3*pow(ax,3));
  
  qanal[anal].qB[bin][j]->Add(qxt, pow(bx,3) - 3*bx*pow(by,2));
  qanal[anal].qB[bin][j]->Add(qyt, pow(by,3) - 3*by*pow(bx,2));
  qanal[anal].qB[bin][j]->Add(qxxy, 9*pow(bx,2)*by - 3*pow(by,3));
  qanal[anal].qB[bin][j]->Add(qxyy, 9*pow(by,2)*bx - 3*pow(bx,3));
  
  qanal[anal].wnA[bin][j]->Add(qct,pow(wa,3));
  qanal[anal].wnB[bin][j]->Add(qct,pow(wb,3));
  qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,3)*pow(bx,3) - 2*pow(ax,2)*ay*pow(bx,3) - ax*pow(ay,2)*pow(bx,3) + 3*pow(ax,3)*pow(bx,2)*by + 6*pow(ax,2)*ay*pow(bx,2)*by - 
		    3*ax*pow(ay,2)*pow(bx,2)*by - 3*pow(ax,3)*bx*pow(by,2) + 6*pow(ax,2)*ay*bx*pow(by,2) + 3*ax*pow(ay,2)*bx*pow(by,2) - 
		    pow(ax,3)*pow(by,3) - 2*pow(ax,2)*ay*pow(by,3) + ax*pow(ay,2)*pow(by,3));
  qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,3)*pow(cx,3) - 2*pow(ax,2)*ay*pow(cx,3) - ax*pow(ay,2)*pow(cx,3) + 3*pow(ax,3)*pow(cx,2)*cy + 6*pow(ax,2)*ay*pow(cx,2)*cy - 
		    3*ax*pow(ay,2)*pow(cx,2)*cy - 3*pow(ax,3)*cx*pow(cy,2) + 6*pow(ax,2)*ay*cx*pow(cy,2) + 3*ax*pow(ay,2)*cx*pow(cy,2) - 
		    pow(ax,3)*pow(cy,3) - 2*pow(ax,2)*ay*pow(cy,3) + ax*pow(ay,2)*pow(cy,3));
  qanal[anal].qCB[bin][j]->Fill(0.,pow(bx,3)*pow(cx,3) - 2*pow(bx,2)*by*pow(cx,3) - bx*pow(by,2)*pow(cx,3) + 3*pow(bx,3)*pow(cx,2)*cy + 6*pow(bx,2)*by*pow(cx,2)*cy - 
		    3*bx*pow(by,2)*pow(cx,2)*cy - 3*pow(bx,3)*cx*pow(cy,2) + 6*pow(bx,2)*by*cx*pow(cy,2) + 3*bx*pow(by,2)*cx*pow(cy,2) - 
		    pow(bx,3)*pow(cy,3) - 2*pow(bx,2)*by*pow(cy,3) + bx*pow(by,2)*pow(cy,3));
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wa*wb,3));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wa*wc,3));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wb*wc,3));
  
}


void Fill_D26A(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp2];
  double ay = qy[HFp2];
  double wa = sumw[HFp2];

  double bx = qx[HFm2];
  double by = qy[HFm2];
  double wb = sumw[HFm2];

  double cx = qx[trackmid2];
  double cy = qy[trackmid2];
  double wc = sumw[trackmid2];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;

  qxxy->Reset();
  qxxy->Add(qxtrk_,1.);
  qxxy->Multiply(qxtrk_);
  qxxy->Multiply(qytrk_);

  qxyy->Reset();
  qxyy->Add(qxtrk_,1.);
  qxyy->Multiply(qytrk_);
  qxyy->Multiply(qytrk_);

  qxt->Reset();
  qxt->Add(qxtrk_,1.);
  qxt->Multiply(qxtrk_);
  qxt->Multiply(qxtrk_);
  qyt->Reset();
  qyt->Add(qytrk_,1.);
  qyt->Multiply(qytrk_);
  qyt->Multiply(qytrk_);
  qct->Reset();
  qct->Add(qcnt_,1.);
  qct->Multiply(qcnt_);
  qct->Multiply(qcnt_);


  qanal[anal].qA[bin][0]->Add(qxt, pow(ax,2)*bx - pow(ay,2)*bx - 2*ax*ay*by);
  qanal[anal].qA[bin][0]->Add(qyt, - 2*ax*ay*bx - pow(ax,2)*by + pow(ay,2)*by);
  qanal[anal].qA[bin][0]->Add(qxxy, 6*ax*ay*bx +  3*pow(ax,2)*by - 3*pow(ay,2)*by);
  qanal[anal].qA[bin][0]->Add(qxyy, - 3*pow(ax,2)*bx +  3*pow(ay,2)*bx + 6*ax*ay*by);
  
  qanal[anal].qB[bin][0]->Add(qxt, pow(bx,2)*ax - pow(by,2)*ax - 2*bx*by*ay);
  qanal[anal].qB[bin][0]->Add(qyt, - 2*bx*by*ax - pow(bx,2)*ay + pow(by,2)*ay);
  qanal[anal].qB[bin][0]->Add(qxxy, 6*bx*by*ax +  3*pow(bx,2)*ay - 3*pow(by,2)*ay);
  qanal[anal].qB[bin][0]->Add(qxyy, - 3*pow(bx,2)*ax +  3*pow(by,2)*ax + 6*bx*by*ay);
  
  qanal[anal].wnA[bin][0]->Add(qct,pow(wa,3));
  qanal[anal].wnB[bin][0]->Add(qct,pow(wb,3));
  qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,3)*pow(bx,3) - 2*pow(ax,2)*ay*pow(bx,3) - ax*pow(ay,2)*pow(bx,3) + 3*pow(ax,3)*pow(bx,2)*by + 6*pow(ax,2)*ay*pow(bx,2)*by - 
		    3*ax*pow(ay,2)*pow(bx,2)*by - 3*pow(ax,3)*bx*pow(by,2) + 6*pow(ax,2)*ay*bx*pow(by,2) + 3*ax*pow(ay,2)*bx*pow(by,2) - 
		    pow(ax,3)*pow(by,3) - 2*pow(ax,2)*ay*pow(by,3) + ax*pow(ay,2)*pow(by,3));
  qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,3)*pow(cx,3) - 2*pow(ax,2)*ay*pow(cx,3) - ax*pow(ay,2)*pow(cx,3) + 3*pow(ax,3)*pow(cx,2)*cy + 6*pow(ax,2)*ay*pow(cx,2)*cy - 
		    3*ax*pow(ay,2)*pow(cx,2)*cy - 3*pow(ax,3)*cx*pow(cy,2) + 6*pow(ax,2)*ay*cx*pow(cy,2) + 3*ax*pow(ay,2)*cx*pow(cy,2) - 
		    pow(ax,3)*pow(cy,3) - 2*pow(ax,2)*ay*pow(cy,3) + ax*pow(ay,2)*pow(cy,3));
  qanal[anal].qCB[bin][0]->Fill(0.,pow(bx,3)*pow(cx,3) - 2*pow(bx,2)*by*pow(cx,3) - bx*pow(by,2)*pow(cx,3) + 3*pow(bx,3)*pow(cx,2)*cy + 6*pow(bx,2)*by*pow(cx,2)*cy - 
		    3*bx*pow(by,2)*pow(cx,2)*cy - 3*pow(bx,3)*cx*pow(cy,2) + 6*pow(bx,2)*by*cx*pow(cy,2) + 3*bx*pow(by,2)*cx*pow(cy,2) - 
		    pow(bx,3)*pow(cy,3) - 2*pow(bx,2)*by*pow(cy,3) + bx*pow(by,2)*pow(cy,3));
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wa*wb,3));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wa*wc,3));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wb*wc,3));
  
  int j=(int)(ran->Uniform(0,9.999))+1;


  qanal[anal].qA[bin][j]->Add(qxt, pow(ax,2)*bx - pow(ay,2)*bx - 2*ax*ay*by);
  qanal[anal].qA[bin][j]->Add(qyt, - 2*ax*ay*bx - pow(ax,2)*by + pow(ay,2)*by);
  qanal[anal].qA[bin][j]->Add(qxxy, 6*ax*ay*bx +  3*pow(ax,2)*by - 3*pow(ay,2)*by);
  qanal[anal].qA[bin][j]->Add(qxyy, - 3*pow(ax,2)*bx +  3*pow(ay,2)*bx + 6*ax*ay*by);
  
  qanal[anal].qB[bin][j]->Add(qxt, pow(bx,2)*ax - pow(by,2)*ax - 2*bx*by*ay);
  qanal[anal].qB[bin][j]->Add(qyt, - 2*bx*by*ax - pow(bx,2)*ay + pow(by,2)*ay);
  qanal[anal].qB[bin][j]->Add(qxxy, 6*bx*by*ax +  3*pow(bx,2)*ay - 3*pow(by,2)*ay);
  qanal[anal].qB[bin][j]->Add(qxyy, - 3*pow(bx,2)*ax +  3*pow(by,2)*ax + 6*bx*by*ay);
  
  qanal[anal].wnA[bin][j]->Add(qct,pow(wa,3));
  qanal[anal].wnB[bin][j]->Add(qct,pow(wb,3));
  qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,3)*pow(bx,3) - 2*pow(ax,2)*ay*pow(bx,3) - ax*pow(ay,2)*pow(bx,3) + 3*pow(ax,3)*pow(bx,2)*by + 6*pow(ax,2)*ay*pow(bx,2)*by - 
		    3*ax*pow(ay,2)*pow(bx,2)*by - 3*pow(ax,3)*bx*pow(by,2) + 6*pow(ax,2)*ay*bx*pow(by,2) + 3*ax*pow(ay,2)*bx*pow(by,2) - 
		    pow(ax,3)*pow(by,3) - 2*pow(ax,2)*ay*pow(by,3) + ax*pow(ay,2)*pow(by,3));
  qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,3)*pow(cx,3) - 2*pow(ax,2)*ay*pow(cx,3) - ax*pow(ay,2)*pow(cx,3) + 3*pow(ax,3)*pow(cx,2)*cy + 6*pow(ax,2)*ay*pow(cx,2)*cy - 
		    3*ax*pow(ay,2)*pow(cx,2)*cy - 3*pow(ax,3)*cx*pow(cy,2) + 6*pow(ax,2)*ay*cx*pow(cy,2) + 3*ax*pow(ay,2)*cx*pow(cy,2) - 
		    pow(ax,3)*pow(cy,3) - 2*pow(ax,2)*ay*pow(cy,3) + ax*pow(ay,2)*pow(cy,3));
  qanal[anal].qCB[bin][j]->Fill(0.,pow(bx,3)*pow(cx,3) - 2*pow(bx,2)*by*pow(cx,3) - bx*pow(by,2)*pow(cx,3) + 3*pow(bx,3)*pow(cx,2)*cy + 6*pow(bx,2)*by*pow(cx,2)*cy - 
		    3*bx*pow(by,2)*pow(cx,2)*cy - 3*pow(bx,3)*cx*pow(cy,2) + 6*pow(bx,2)*by*cx*pow(cy,2) + 3*bx*pow(by,2)*cx*pow(cy,2) - 
		    pow(bx,3)*pow(cy,3) - 2*pow(bx,2)*by*pow(cy,3) + bx*pow(by,2)*pow(cy,3));
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wa*wb,3));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wa*wc,3));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wb*wc,3));
  
}
void Fill_D34(int anal, int bin, TH2D * qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double ax = qx[HFp3];
  double ay = qy[HFp3];
  double wa = sumw[HFp3];

  double bx = qx[HFm3];
  double by = qy[HFm3];
  double wb = sumw[HFm3];

  double cx = qx[trackmid3];
  double cy = qy[trackmid3];
  double wc = sumw[trackmid3];

  if(pow(ax,2)+pow(ay,2) < 1e-6) return;
  if(pow(bx,2)+pow(by,2) < 1e-6) return;
  if(pow(cx,2)+pow(cy,2) < 1e-6) return;

  qxy->Reset();
  qxy->Add(qxtrk_,1.);
  qxy->Multiply(qytrk_);

  qxt->Reset();
  qxt->Add(qxtrk_,1.);
  qxt->Multiply(qxtrk_);
  qyt->Reset();
  qyt->Add(qytrk_,1.);
  qyt->Multiply(qytrk_);
  qct->Reset();
  qct->Add(qcnt_,1.);
  qct->Multiply(qcnt_);


  qanal[anal].qA[bin][0]->Add(qxt, pow(ax,2) - pow(ay,2));
  qanal[anal].qA[bin][0]->Add(qyt, pow(ay,2) - pow(ax,2));
  qanal[anal].qA[bin][0]->Add(qxy, 4*ax*ay);

  qanal[anal].qB[bin][0]->Add(qxt, pow(bx,2) - pow(by,2));
  qanal[anal].qB[bin][0]->Add(qyt, pow(by,2) - pow(bx,2));
  qanal[anal].qB[bin][0]->Add(qxy, 4*bx*by);

  qanal[anal].wnA[bin][0]->Add(qct,pow(wa,2));
  qanal[anal].wnB[bin][0]->Add(qct,pow(wb,2));
  qanal[anal].qBA[bin][0]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,pow(wa*wb,2));
  qanal[anal].qCAcnt[bin][0]->Fill(0.,pow(wa*wc,2));
  qanal[anal].qCBcnt[bin][0]->Fill(0.,pow(wb*wc,2));

  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxt, pow(ax,2) - pow(ay,2));
  qanal[anal].qA[bin][j]->Add(qyt, pow(ay,2) - pow(ax,2));
  qanal[anal].qA[bin][j]->Add(qxy, 4*ax*ay);
  
  qanal[anal].qB[bin][j]->Add(qxt, pow(bx,2) - pow(by,2));
  qanal[anal].qB[bin][j]->Add(qyt, pow(by,2) - pow(bx,2));
  qanal[anal].qB[bin][j]->Add(qxy, 4*bx*by);
  
  qanal[anal].wnA[bin][j]->Add(qct,pow(wa,2));
  qanal[anal].wnB[bin][j]->Add(qct,pow(wb,2));
  qanal[anal].qBA[bin][j]->Fill(0.,pow(ax,2)*(pow(bx,2)-pow(by,2)) + pow(ay,2)*(pow(by,2)-pow(bx,2))+4*ax*ay*bx*by);
  qanal[anal].qCA[bin][j]->Fill(0.,pow(ax,2)*(pow(cx,2)-pow(cy,2)) + pow(ay,2)*(pow(cy,2)-pow(cx,2))+4*ax*ay*cx*cy);
  qanal[anal].qCB[bin][j]->Fill(0.,pow(bx,2)*(pow(cx,2)-pow(cy,2)) + pow(by,2)*(pow(cy,2)-pow(cx,2))+4*bx*by*cx*cy);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,pow(wa*wb,2));
  qanal[anal].qCAcnt[bin][j]->Fill(0.,pow(wa*wc,2));
  qanal[anal].qCBcnt[bin][j]->Fill(0.,pow(wb*wc,2));
  
  
}

void Fill_D34A(int anal, int bin, TH2D * qxtrk2_, TH2D * qytrk2_, TH2D * qxtrk3_, TH2D * qytrk3_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp3];
  double a2y = qy[HFp3];
  double a2w = sumw[HFp3];
  double a3x = qx[HFm3];
  double a3y = qy[HFm3];
  double a3w = sumw[HFm3];
  
  double b2x = qx[HFm3];
  double b2y = qy[HFm3];
  double b2w = sumw[HFm3];
  double b3x = qx[HFp3];
  double b3y = qy[HFp3];
  double b3w = sumw[HFp3];
  
  double c2x = qx[trackmid3];
  double c2y = qy[trackmid3];
  double c2w = sumw[trackmid3];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];

  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  qx2y3->Reset();
  qx2y3->Add(qxtrk2_,1.);
  qx2y3->Multiply(qytrk3_);

  qy2x3->Reset();
  qy2x3->Add(qytrk2_,1.);
  qy2x3->Multiply(qxtrk3_);


  qxt->Reset();
  qxt->Add(qxtrk2_,1.);
  qxt->Multiply(qxtrk3_);
  qyt->Reset();
  qyt->Add(qytrk2_,1.);
  qyt->Multiply(qytrk3_);
  qct->Reset();
  qct->Add(qcnt_,1.);
  qct->Multiply(qcnt_);

  qanal[anal].qA[bin][0]->Add(qxt,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][0]->Add(qyt,a2y*a3y-a2x*a3x);
  qanal[anal].qA[bin][0]->Add(qy2x3,  a2y*a3x+a2x*a3y);
  qanal[anal].qA[bin][0]->Add(qx2y3,  a2y*a3x+a2x*a3y);

  qanal[anal].qB[bin][0]->Add(qxt,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][0]->Add(qyt,b2y*b3y-b2x*b3x);
  qanal[anal].qB[bin][0]->Add(qy2x3,  b2y*b3x+b2x*b3y);
  qanal[anal].qB[bin][0]->Add(qx2y3,  b2y*b3x+b2x*b3y);

  qanal[anal].wnA[bin][0]->Add(qct,a2w*a3w);
  qanal[anal].wnB[bin][0]->Add(qct,b2w*b3w);

  qanal[anal].qBA[bin][0]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c3w*b2w*b3w);
  
  int j=(int)(ran->Uniform(0,9.999))+1;

  qanal[anal].qA[bin][j]->Add(qxt,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][j]->Add(qyt,a2y*a3y-a2x*a3x);
  qanal[anal].qA[bin][j]->Add(qy2x3,  a2y*a3x+a2x*a3y);
  qanal[anal].qA[bin][j]->Add(qx2y3,  a2y*a3x+a2x*a3y);

  qanal[anal].qB[bin][j]->Add(qxt,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][j]->Add(qyt,b2y*b3y-b2x*b3x);
  qanal[anal].qB[bin][j]->Add(qy2x3,  b2y*b3x+b2x*b3y);
  qanal[anal].qB[bin][j]->Add(qx2y3,  b2y*b3x+b2x*b3y);

  qanal[anal].wnA[bin][j]->Add(qct,a2w*a3w);
  qanal[anal].wnB[bin][j]->Add(qct,b2w*b3w);

  qanal[anal].qBA[bin][j]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][j]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][j]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,c2w*c3w*b2w*b3w);
  
  
}

void Fill_D2232(int anal, int bin, TH2D * qxtrk2_, TH2D * qytrk2_, TH2D * qxtrk3_, TH2D * qytrk3_, TH2D * qcnt_, double * qx, double * qy, double * sumw){

  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFp3];
  double a3y = qy[HFp3];
  double a3w = sumw[HFp3];

  double b2x = qx[HFm3];
  double b2y = qy[HFm3];
  double b2w = sumw[HFm3];
  double b3x = qx[HFm3];
  double b3y = qy[HFm3];
  double b3w = sumw[HFm3];

  double c2x = qx[trackmid3];
  double c2y = qy[trackmid3];
  double c2w = sumw[trackmid3];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];

  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  qx2y3->Reset();
  qx2y3->Add(qxtrk2_,1.);
  qx2y3->Multiply(qytrk3_);

  qy2x3->Reset();
  qy2x3->Add(qytrk2_,1.);
  qy2x3->Multiply(qxtrk3_);

  qxt->Reset();
  qxt->Add(qxtrk2_,1.);
  qxt->Multiply(qxtrk3_);
  qyt->Reset();
  qyt->Add(qytrk2_,1.);
  qyt->Multiply(qytrk3_);
  qct->Reset();
  qct->Add(qcnt_,1.);
  qct->Multiply(qcnt_);


  //  qxtrk2_->Multiply(qxtrk3_);
  //  qytrk2_->Multiply(qytrk3_);
  //
  //qcnt_->Multiply(qcnt_);

  qanal[anal].qA[bin][0]->Add(qxt,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][0]->Add(qyt,a2y*a3y-a2x*a3x);
  qanal[anal].qA[bin][0]->Add(qy2x3,  a2y*a3x+a2x*a3y);
  qanal[anal].qA[bin][0]->Add(qx2y3,  a2y*a3x+a2x*a3y);

  qanal[anal].qB[bin][0]->Add(qxt,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][0]->Add(qyt,b2y*b3y-b2x*b3x);
  qanal[anal].qB[bin][0]->Add(qy2x3,  b2y*b3x+b2x*b3y);
  qanal[anal].qB[bin][0]->Add(qx2y3,  b2y*b3x+b2x*b3y);

  qanal[anal].wnA[bin][0]->Add(qct,a2w*a3w);
  qanal[anal].wnB[bin][0]->Add(qct,b2w*b3w);

  qanal[anal].qBA[bin][0]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c3w*b2w*b3w);
  
  int j=(int)(ran->Uniform(0,9.999))+1;


  qanal[anal].qA[bin][j]->Add(qxt,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][j]->Add(qyt,a2y*a3y-a2x*a3x);
  qanal[anal].qA[bin][j]->Add(qy2x3,  a2y*a3x+a2x*a3y);
  qanal[anal].qA[bin][j]->Add(qx2y3,  a2y*a3x+a2x*a3y);

  qanal[anal].qB[bin][j]->Add(qxt,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][j]->Add(qyt,b2y*b3y-b2x*b3x);
  qanal[anal].qB[bin][j]->Add(qy2x3,  b2y*b3x+b2x*b3y);
  qanal[anal].qB[bin][j]->Add(qx2y3,  b2y*b3x+b2x*b3y);

  qanal[anal].wnA[bin][j]->Add(qct,a2w*a3w);
  qanal[anal].wnB[bin][j]->Add(qct,b2w*b3w);

  qanal[anal].qBA[bin][j]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][j]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][j]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,c2w*c3w*b2w*b3w);
  
  
}

void Fill_D2232A(int anal, int bin, TH2D * qxtrk2_, TH2D * qytrk2_, TH2D * qxtrk3_, TH2D * qytrk3_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFm3];
  double a3y = qy[HFm3];
  double a3w = sumw[HFm3];
  
  double b2x = qx[HFm2];
  double b2y = qy[HFm2];
  double b2w = sumw[HFm2];
  double b3x = qx[HFp3];
  double b3y = qy[HFp3];
  double b3w = sumw[HFp3];
  
  double c2x = qx[trackmid2];
  double c2y = qy[trackmid2];
  double c2w = sumw[trackmid2];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];


  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  qx2y3->Reset();
  qx2y3->Add(qxtrk2_,1.);
  qx2y3->Multiply(qytrk3_);

  qy2x3->Reset();
  qy2x3->Add(qytrk2_,1.);
  qy2x3->Multiply(qxtrk3_);


  qxt->Reset();
  qxt->Add(qxtrk2_,1.);
  qxt->Multiply(qxtrk3_);
  qyt->Reset();
  qyt->Add(qytrk2_,1.);
  qyt->Multiply(qytrk3_);
  qct->Reset();
  qct->Add(qcnt_,1.);
  qct->Multiply(qcnt_);

  qanal[anal].qA[bin][0]->Add(qxt,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][0]->Add(qyt,a2y*a3y-a2x*a3x);
  qanal[anal].qA[bin][0]->Add(qy2x3,  a2y*a3x+a2x*a3y);
  qanal[anal].qA[bin][0]->Add(qx2y3,  a2y*a3x+a2x*a3y);

  qanal[anal].qB[bin][0]->Add(qxt,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][0]->Add(qyt,b2y*b3y-b2x*b3x);
  qanal[anal].qB[bin][0]->Add(qy2x3,  b2y*b3x+b2x*b3y);
  qanal[anal].qB[bin][0]->Add(qx2y3,  b2y*b3x+b2x*b3y);

  qanal[anal].wnA[bin][0]->Add(qct,a2w*a3w);
  qanal[anal].wnB[bin][0]->Add(qct,b2w*b3w);

  qanal[anal].qBA[bin][0]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c3w*b2w*b3w);
  
  int j=(int)(ran->Uniform(0,9.999))+1;

  qanal[anal].qA[bin][j]->Add(qxt,a2x*a3x-a2y*a3y);
  qanal[anal].qA[bin][j]->Add(qyt,a2y*a3y-a2x*a3x);
  qanal[anal].qA[bin][j]->Add(qy2x3,  a2y*a3x+a2x*a3y);
  qanal[anal].qA[bin][j]->Add(qx2y3,  a2y*a3x+a2x*a3y);

  qanal[anal].qB[bin][j]->Add(qxt,b2x*b3x-b2y*b3y);
  qanal[anal].qB[bin][j]->Add(qyt,b2y*b3y-b2x*b3x);
  qanal[anal].qB[bin][j]->Add(qy2x3,  b2y*b3x+b2x*b3y);
  qanal[anal].qB[bin][j]->Add(qx2y3,  b2y*b3x+b2x*b3y);

  qanal[anal].wnA[bin][j]->Add(qct,a2w*a3w);
  qanal[anal].wnB[bin][j]->Add(qct,b2w*b3w);

  qanal[anal].qBA[bin][j]->Fill(0.,a2x*a3x*b2x*b3x - a2y*a3y*b2x*b3x + a2y*a3x*b2y*b3x + a2x*a3y*b2y*b3x + a2y*a3x*b2x*b3y + a2x*a3y*b2x*b3y - a2x*a3x*b2y*b3y + a2y*a3y*b2y*b3y);
  qanal[anal].qCA[bin][j]->Fill(0.,a2x*a3x*c2x*c3x - a2y*a3y*c2x*c3x + a2y*a3x*c2y*c3x + a2x*a3y*c2y*c3x + a2y*a3x*c2x*c3y + a2x*a3y*c2x*c3y - a2x*a3x*c2y*c3y + a2y*a3y*c2y*c3y);
  qanal[anal].qCB[bin][j]->Fill(0.,b2x*b3x*c2x*c3x - b2y*b3y*c2x*c3x + b2y*b3x*c2y*c3x + b2x*b3y*c2y*c3x + b2y*b3x*c2x*c3y + b2x*b3y*c2x*c3y - b2x*b3x*c2y*c3y + b2y*b3y*c2y*c3y);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,b2w*b3w*a2w*a3w);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,c2w*c3w*a2w*a3w);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,c2w*c3w*b2w*b3w);
  
  
}


void Fill_D2432(int anal, int bin, TH2D * qxtrk2_, TH2D * qytrk2_, TH2D * qxtrk3_, TH2D * qytrk3_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFp3];
  double a3y = qy[HFp3];
  double a3w = sumw[HFp3];
  
  double b2x = qx[HFm2];
  double b2y = qy[HFm2];
  double b2w = sumw[HFm2];
  double b3x = qx[HFm3];
  double b3y = qy[HFm3];
  double b3w = sumw[HFm3];
  
  double c2x = qx[trackmid2];
  double c2y = qy[trackmid2];
  double c2w = sumw[trackmid2];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];

  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  int k=(int)(ran->Uniform(0,9.999))+1;
  for(int i = 1; i<=qxtrk2_->GetNbinsX(); i++) {
    for(int j = 1; j<=qxtrk2_->GetNbinsY(); j++) {
      double q2x = qxtrk2_->GetBinContent(i,j);
      double q2y = qytrk2_->GetBinContent(i,j);
      double q3x = qxtrk3_->GetBinContent(i,j);
      double q3y = qytrk3_->GetBinContent(i,j);
      double cnt = qcnt_->GetBinContent(i,j);

      double qAterm = pow(a2x,2)*a3x*pow(q2x,2)*q3x - pow(a2y,2)*a3x*pow(q2x,2)*q3x - 2*a2x*a2y*a3y*pow(q2x,2)*q3x + 4*a2x*a2y*a3x*q2x*q2y*q3x + 
	2*pow(a2x,2)*a3y*q2x*q2y*q3x - 2*pow(a2y,2)*a3y*q2x*q2y*q3x - pow(a2x,2)*a3x*pow(q2y,2)*q3x + pow(a2y,2)*a3x*pow(q2y,2)*q3x + 
	2*a2x*a2y*a3y*pow(q2y,2)*q3x + 2*a2x*a2y*a3x*pow(q2x,2)*q3y + pow(a2x,2)*a3y*pow(q2x,2)*q3y - pow(a2y,2)*a3y*pow(q2x,2)*q3y - 
	2*pow(a2x,2)*a3x*q2x*q2y*q3y + 2*pow(a2y,2)*a3x*q2x*q2y*q3y + 4*a2x*a2y*a3y*q2x*q2y*q3y - 2*a2x*a2y*a3x*pow(q2y,2)*q3y - 
	pow(a2x,2)*a3y*pow(q2y,2)*q3y + pow(a2y,2)*a3y*pow(q2y,2)*q3y;
      double qBterm = pow(b2x,2)*b3x*pow(q2x,2)*q3x - pow(b2y,2)*b3x*pow(q2x,2)*q3x - 2*b2x*b2y*b3y*pow(q2x,2)*q3x + 4*b2x*b2y*b3x*q2x*q2y*q3x + 
	2*pow(b2x,2)*b3y*q2x*q2y*q3x - 2*pow(b2y,2)*b3y*q2x*q2y*q3x - pow(b2x,2)*b3x*pow(q2y,2)*q3x + pow(b2y,2)*b3x*pow(q2y,2)*q3x + 
	2*b2x*b2y*b3y*pow(q2y,2)*q3x + 2*b2x*b2y*b3x*pow(q2x,2)*q3y + pow(b2x,2)*b3y*pow(q2x,2)*q3y - pow(b2y,2)*b3y*pow(q2x,2)*q3y - 
	2*pow(b2x,2)*b3x*q2x*q2y*q3y + 2*pow(b2y,2)*b3x*q2x*q2y*q3y + 4*b2x*b2y*b3y*q2x*q2y*q3y - 2*b2x*b2y*b3x*pow(q2y,2)*q3y - 
	pow(b2x,2)*b3y*pow(q2y,2)*q3y + pow(b2y,2)*b3y*pow(q2y,2)*q3y;
      qanal[anal].qA[bin][0]->SetBinContent(i,j,qanal[anal].qA[bin][0]->GetBinContent(i,j) + qAterm);
      qanal[anal].qB[bin][0]->SetBinContent(i,j,qanal[anal].qB[bin][0]->GetBinContent(i,j) + qBterm);
      qanal[anal].wnA[bin][0]->SetBinContent(i,j,qanal[anal].wnA[bin][0]->GetBinContent(i,j) + pow(cnt,3)*a2w*a2w*a3w);
      qanal[anal].wnB[bin][0]->SetBinContent(i,j,qanal[anal].wnB[bin][0]->GetBinContent(i,j) + pow(cnt,3)*b2w*b2w*b3w);
      qanal[anal].qA[bin][k]->SetBinContent(i,j,qanal[anal].qA[bin][k]->GetBinContent(i,j) + qAterm);
      qanal[anal].qB[bin][k]->SetBinContent(i,j,qanal[anal].qB[bin][k]->GetBinContent(i,j) + qBterm);
      qanal[anal].wnA[bin][k]->SetBinContent(i,j,qanal[anal].wnA[bin][k]->GetBinContent(i,j) + pow(cnt,3)*a2w*a2w*a3w);
      qanal[anal].wnB[bin][k]->SetBinContent(i,j,qanal[anal].wnB[bin][k]->GetBinContent(i,j) + pow(cnt,3)*b2w*b2w*b3w);
      
    }
  }

  qanal[anal].qBA[bin][0]->Fill(0.,pow(a2x,2)*a3x*pow(b2x,2)*b3x - pow(a2y,2)*a3x*pow(b2x,2)*b3x - 2*a2x*a2y*a3y*pow(b2x,2)*b3x + 4*a2x*a2y*a3x*b2x*b2y*b3x + 2*pow(a2x,2)*a3y*b2x*b2y*b3x - 
		    2*pow(a2y,2)*a3y*b2x*b2y*b3x - pow(a2x,2)*a3x*pow(b2y,2)*b3x + pow(a2y,2)*a3x*pow(b2y,2)*b3x + 2*a2x*a2y*a3y*pow(b2y,2)*b3x + 2*a2x*a2y*a3x*pow(b2x,2)*b3y + 
		    pow(a2x,2)*a3y*pow(b2x,2)*b3y - pow(a2y,2)*a3y*pow(b2x,2)*b3y - 2*pow(a2x,2)*a3x*b2x*b2y*b3y + 2*pow(a2y,2)*a3x*b2x*b2y*b3y + 4*a2x*a2y*a3y*b2x*b2y*b3y - 
		    2*a2x*a2y*a3x*pow(b2y,2)*b3y - pow(a2x,2)*a3y*pow(b2y,2)*b3y + pow(a2y,2)*a3y*pow(b2y,2)*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(a2x,2)*a3x*pow(c2x,2)*c3x - pow(a2y,2)*a3x*pow(c2x,2)*c3x - 2*a2x*a2y*a3y*pow(c2x,2)*c3x + 4*a2x*a2y*a3x*c2x*c2y*c3x + 2*pow(a2x,2)*a3y*c2x*c2y*c3x - 
		    2*pow(a2y,2)*a3y*c2x*c2y*c3x - pow(a2x,2)*a3x*pow(c2y,2)*c3x + pow(a2y,2)*a3x*pow(c2y,2)*c3x + 2*a2x*a2y*a3y*pow(c2y,2)*c3x + 2*a2x*a2y*a3x*pow(c2x,2)*c3y + 
		    pow(a2x,2)*a3y*pow(c2x,2)*c3y - pow(a2y,2)*a3y*pow(c2x,2)*c3y - 2*pow(a2x,2)*a3x*c2x*c2y*c3y + 2*pow(a2y,2)*a3x*c2x*c2y*c3y + 4*a2x*a2y*a3y*c2x*c2y*c3y - 
		    2*a2x*a2y*a3x*pow(c2y,2)*c3y - pow(a2x,2)*a3y*pow(c2y,2)*c3y + pow(a2y,2)*a3y*pow(c2y,2)*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(b2x,2)*b3x*pow(c2x,2)*c3x - pow(b2y,2)*b3x*pow(c2x,2)*c3x - 2*b2x*b2y*b3y*pow(c2x,2)*c3x + 4*b2x*b2y*b3x*c2x*c2y*c3x + 2*pow(b2x,2)*b3y*c2x*c2y*c3x - 
		    2*pow(b2y,2)*b3y*c2x*c2y*c3x - pow(b2x,2)*b3x*pow(c2y,2)*c3x + pow(b2y,2)*b3x*pow(c2y,2)*c3x + 2*b2x*b2y*b3y*pow(c2y,2)*c3x + 2*b2x*b2y*b3x*pow(c2x,2)*c3y + 
		    pow(b2x,2)*b3y*pow(c2x,2)*c3y - pow(b2y,2)*b3y*pow(c2x,2)*c3y - 2*pow(b2x,2)*b3x*c2x*c2y*c3y + 2*pow(b2y,2)*b3x*c2x*c2y*c3y + 4*b2x*b2y*b3y*c2x*c2y*c3y - 
		    2*b2x*b2y*b3x*pow(c2y,2)*c3y - pow(b2x,2)*b3y*pow(c2y,2)*c3y + pow(b2y,2)*b3y*pow(c2y,2)*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b2w*b3w*a2w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c2w*c3w*a2w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c2w*c3w*b2w*b2w*b3w);

  qanal[anal].qBA[bin][k]->Fill(0.,pow(a2x,2)*a3x*pow(b2x,2)*b3x - pow(a2y,2)*a3x*pow(b2x,2)*b3x - 2*a2x*a2y*a3y*pow(b2x,2)*b3x + 4*a2x*a2y*a3x*b2x*b2y*b3x + 2*pow(a2x,2)*a3y*b2x*b2y*b3x - 
		    2*pow(a2y,2)*a3y*b2x*b2y*b3x - pow(a2x,2)*a3x*pow(b2y,2)*b3x + pow(a2y,2)*a3x*pow(b2y,2)*b3x + 2*a2x*a2y*a3y*pow(b2y,2)*b3x + 2*a2x*a2y*a3x*pow(b2x,2)*b3y + 
		    pow(a2x,2)*a3y*pow(b2x,2)*b3y - pow(a2y,2)*a3y*pow(b2x,2)*b3y - 2*pow(a2x,2)*a3x*b2x*b2y*b3y + 2*pow(a2y,2)*a3x*b2x*b2y*b3y + 4*a2x*a2y*a3y*b2x*b2y*b3y - 
		    2*a2x*a2y*a3x*pow(b2y,2)*b3y - pow(a2x,2)*a3y*pow(b2y,2)*b3y + pow(a2y,2)*a3y*pow(b2y,2)*b3y);
  qanal[anal].qCA[bin][k]->Fill(0.,pow(a2x,2)*a3x*pow(c2x,2)*c3x - pow(a2y,2)*a3x*pow(c2x,2)*c3x - 2*a2x*a2y*a3y*pow(c2x,2)*c3x + 4*a2x*a2y*a3x*c2x*c2y*c3x + 2*pow(a2x,2)*a3y*c2x*c2y*c3x - 
		    2*pow(a2y,2)*a3y*c2x*c2y*c3x - pow(a2x,2)*a3x*pow(c2y,2)*c3x + pow(a2y,2)*a3x*pow(c2y,2)*c3x + 2*a2x*a2y*a3y*pow(c2y,2)*c3x + 2*a2x*a2y*a3x*pow(c2x,2)*c3y + 
		    pow(a2x,2)*a3y*pow(c2x,2)*c3y - pow(a2y,2)*a3y*pow(c2x,2)*c3y - 2*pow(a2x,2)*a3x*c2x*c2y*c3y + 2*pow(a2y,2)*a3x*c2x*c2y*c3y + 4*a2x*a2y*a3y*c2x*c2y*c3y - 
		    2*a2x*a2y*a3x*pow(c2y,2)*c3y - pow(a2x,2)*a3y*pow(c2y,2)*c3y + pow(a2y,2)*a3y*pow(c2y,2)*c3y);
  qanal[anal].qCB[bin][k]->Fill(0.,pow(b2x,2)*b3x*pow(c2x,2)*c3x - pow(b2y,2)*b3x*pow(c2x,2)*c3x - 2*b2x*b2y*b3y*pow(c2x,2)*c3x + 4*b2x*b2y*b3x*c2x*c2y*c3x + 2*pow(b2x,2)*b3y*c2x*c2y*c3x - 
		    2*pow(b2y,2)*b3y*c2x*c2y*c3x - pow(b2x,2)*b3x*pow(c2y,2)*c3x + pow(b2y,2)*b3x*pow(c2y,2)*c3x + 2*b2x*b2y*b3y*pow(c2y,2)*c3x + 2*b2x*b2y*b3x*pow(c2x,2)*c3y + 
		    pow(b2x,2)*b3y*pow(c2x,2)*c3y - pow(b2y,2)*b3y*pow(c2x,2)*c3y - 2*pow(b2x,2)*b3x*c2x*c2y*c3y + 2*pow(b2y,2)*b3x*c2x*c2y*c3y + 4*b2x*b2y*b3y*c2x*c2y*c3y - 
		    2*b2x*b2y*b3x*pow(c2y,2)*c3y - pow(b2x,2)*b3y*pow(c2y,2)*c3y + pow(b2y,2)*b3y*pow(c2y,2)*c3y);
  qanal[anal].qBAcnt[bin][k]->Fill(0.,b2w*b2w*b3w*a2w*a2w*a3w);
  qanal[anal].qCAcnt[bin][k]->Fill(0.,c2w*c2w*c3w*a2w*a2w*a3w);
  qanal[anal].qCBcnt[bin][k]->Fill(0.,c2w*c2w*c3w*b2w*b2w*b3w);


}

void Fill_D2432A(int anal, int bin, TH2D * qxtrk2_, TH2D * qytrk2_, TH2D * qxtrk3_, TH2D * qytrk3_, TH2D * qcnt_, double * qx, double * qy, double * sumw){
  double a2x = qx[HFp2];
  double a2y = qy[HFp2];
  double a2w = sumw[HFp2];
  double a3x = qx[HFm3];
  double a3y = qy[HFm3];
  double a3w = sumw[HFm3];
  
  double b2x = qx[HFm2];
  double b2y = qy[HFm2];
  double b2w = sumw[HFm2];
  double b3x = qx[HFp3];
  double b3y = qy[HFp3];
  double b3w = sumw[HFp3];
  
  double c2x = qx[trackmid2];
  double c2y = qy[trackmid2];
  double c2w = sumw[trackmid2];
  double c3x = qx[trackmid3];
  double c3y = qy[trackmid3];
  double c3w = sumw[trackmid3];

  if(pow(a2x,2)+pow(a2y,2) < 1e-6) return;
  if(pow(b2x,2)+pow(b2y,2) < 1e-6) return;
  if(pow(c2x,2)+pow(c2y,2) < 1e-6) return;
  if(pow(a3x,2)+pow(a3y,2) < 1e-6) return;
  if(pow(b3x,2)+pow(b3y,2) < 1e-6) return;
  if(pow(c3x,2)+pow(c3y,2) < 1e-6) return;

  int k=(int)(ran->Uniform(0,9.999))+1;
  for(int i = 1; i<=qxtrk2_->GetNbinsX(); i++) {
    for(int j = 1; j<=qxtrk2_->GetNbinsY(); j++) {
      double q2x = qxtrk2_->GetBinContent(i,j);
      double q2y = qytrk2_->GetBinContent(i,j);
      double q3x = qxtrk3_->GetBinContent(i,j);
      double q3y = qytrk3_->GetBinContent(i,j);
      double cnt = qcnt_->GetBinContent(i,j);

      double qAterm = pow(a2x,2)*a3x*pow(q2x,2)*q3x - pow(a2y,2)*a3x*pow(q2x,2)*q3x - 2*a2x*a2y*a3y*pow(q2x,2)*q3x + 4*a2x*a2y*a3x*q2x*q2y*q3x + 
	2*pow(a2x,2)*a3y*q2x*q2y*q3x - 2*pow(a2y,2)*a3y*q2x*q2y*q3x - pow(a2x,2)*a3x*pow(q2y,2)*q3x + pow(a2y,2)*a3x*pow(q2y,2)*q3x + 
	2*a2x*a2y*a3y*pow(q2y,2)*q3x + 2*a2x*a2y*a3x*pow(q2x,2)*q3y + pow(a2x,2)*a3y*pow(q2x,2)*q3y - pow(a2y,2)*a3y*pow(q2x,2)*q3y - 
	2*pow(a2x,2)*a3x*q2x*q2y*q3y + 2*pow(a2y,2)*a3x*q2x*q2y*q3y + 4*a2x*a2y*a3y*q2x*q2y*q3y - 2*a2x*a2y*a3x*pow(q2y,2)*q3y - 
	pow(a2x,2)*a3y*pow(q2y,2)*q3y + pow(a2y,2)*a3y*pow(q2y,2)*q3y;
      double qBterm = pow(b2x,2)*b3x*pow(q2x,2)*q3x - pow(b2y,2)*b3x*pow(q2x,2)*q3x - 2*b2x*b2y*b3y*pow(q2x,2)*q3x + 4*b2x*b2y*b3x*q2x*q2y*q3x + 
	2*pow(b2x,2)*b3y*q2x*q2y*q3x - 2*pow(b2y,2)*b3y*q2x*q2y*q3x - pow(b2x,2)*b3x*pow(q2y,2)*q3x + pow(b2y,2)*b3x*pow(q2y,2)*q3x + 
	2*b2x*b2y*b3y*pow(q2y,2)*q3x + 2*b2x*b2y*b3x*pow(q2x,2)*q3y + pow(b2x,2)*b3y*pow(q2x,2)*q3y - pow(b2y,2)*b3y*pow(q2x,2)*q3y - 
	2*pow(b2x,2)*b3x*q2x*q2y*q3y + 2*pow(b2y,2)*b3x*q2x*q2y*q3y + 4*b2x*b2y*b3y*q2x*q2y*q3y - 2*b2x*b2y*b3x*pow(q2y,2)*q3y - 
	pow(b2x,2)*b3y*pow(q2y,2)*q3y + pow(b2y,2)*b3y*pow(q2y,2)*q3y;
      qanal[anal].qA[bin][0]->SetBinContent(i,j,qanal[anal].qA[bin][0]->GetBinContent(i,j) + qAterm);
      qanal[anal].qB[bin][0]->SetBinContent(i,j,qanal[anal].qB[bin][0]->GetBinContent(i,j) + qBterm);
      qanal[anal].wnA[bin][0]->SetBinContent(i,j,qanal[anal].wnA[bin][0]->GetBinContent(i,j) + pow(cnt,3)*a2w*a2w*a3w);
      qanal[anal].wnB[bin][0]->SetBinContent(i,j,qanal[anal].wnB[bin][0]->GetBinContent(i,j) + pow(cnt,3)*b2w*b2w*b3w);
      qanal[anal].qA[bin][k]->SetBinContent(i,j,qanal[anal].qA[bin][k]->GetBinContent(i,j) + qAterm);
      qanal[anal].qB[bin][k]->SetBinContent(i,j,qanal[anal].qB[bin][k]->GetBinContent(i,j) + qBterm);
      qanal[anal].wnA[bin][k]->SetBinContent(i,j,qanal[anal].wnA[bin][k]->GetBinContent(i,j) + pow(cnt,3)*a2w*a2w*a3w);
      qanal[anal].wnB[bin][k]->SetBinContent(i,j,qanal[anal].wnB[bin][k]->GetBinContent(i,j) + pow(cnt,3)*b2w*b2w*b3w);
      
    }
  }

  qanal[anal].qBA[bin][0]->Fill(0.,pow(a2x,2)*a3x*pow(b2x,2)*b3x - pow(a2y,2)*a3x*pow(b2x,2)*b3x - 2*a2x*a2y*a3y*pow(b2x,2)*b3x + 4*a2x*a2y*a3x*b2x*b2y*b3x + 2*pow(a2x,2)*a3y*b2x*b2y*b3x - 
		    2*pow(a2y,2)*a3y*b2x*b2y*b3x - pow(a2x,2)*a3x*pow(b2y,2)*b3x + pow(a2y,2)*a3x*pow(b2y,2)*b3x + 2*a2x*a2y*a3y*pow(b2y,2)*b3x + 2*a2x*a2y*a3x*pow(b2x,2)*b3y + 
		    pow(a2x,2)*a3y*pow(b2x,2)*b3y - pow(a2y,2)*a3y*pow(b2x,2)*b3y - 2*pow(a2x,2)*a3x*b2x*b2y*b3y + 2*pow(a2y,2)*a3x*b2x*b2y*b3y + 4*a2x*a2y*a3y*b2x*b2y*b3y - 
		    2*a2x*a2y*a3x*pow(b2y,2)*b3y - pow(a2x,2)*a3y*pow(b2y,2)*b3y + pow(a2y,2)*a3y*pow(b2y,2)*b3y);
  qanal[anal].qCA[bin][0]->Fill(0.,pow(a2x,2)*a3x*pow(c2x,2)*c3x - pow(a2y,2)*a3x*pow(c2x,2)*c3x - 2*a2x*a2y*a3y*pow(c2x,2)*c3x + 4*a2x*a2y*a3x*c2x*c2y*c3x + 2*pow(a2x,2)*a3y*c2x*c2y*c3x - 
		    2*pow(a2y,2)*a3y*c2x*c2y*c3x - pow(a2x,2)*a3x*pow(c2y,2)*c3x + pow(a2y,2)*a3x*pow(c2y,2)*c3x + 2*a2x*a2y*a3y*pow(c2y,2)*c3x + 2*a2x*a2y*a3x*pow(c2x,2)*c3y + 
		    pow(a2x,2)*a3y*pow(c2x,2)*c3y - pow(a2y,2)*a3y*pow(c2x,2)*c3y - 2*pow(a2x,2)*a3x*c2x*c2y*c3y + 2*pow(a2y,2)*a3x*c2x*c2y*c3y + 4*a2x*a2y*a3y*c2x*c2y*c3y - 
		    2*a2x*a2y*a3x*pow(c2y,2)*c3y - pow(a2x,2)*a3y*pow(c2y,2)*c3y + pow(a2y,2)*a3y*pow(c2y,2)*c3y);
  qanal[anal].qCB[bin][0]->Fill(0.,pow(b2x,2)*b3x*pow(c2x,2)*c3x - pow(b2y,2)*b3x*pow(c2x,2)*c3x - 2*b2x*b2y*b3y*pow(c2x,2)*c3x + 4*b2x*b2y*b3x*c2x*c2y*c3x + 2*pow(b2x,2)*b3y*c2x*c2y*c3x - 
		    2*pow(b2y,2)*b3y*c2x*c2y*c3x - pow(b2x,2)*b3x*pow(c2y,2)*c3x + pow(b2y,2)*b3x*pow(c2y,2)*c3x + 2*b2x*b2y*b3y*pow(c2y,2)*c3x + 2*b2x*b2y*b3x*pow(c2x,2)*c3y + 
		    pow(b2x,2)*b3y*pow(c2x,2)*c3y - pow(b2y,2)*b3y*pow(c2x,2)*c3y - 2*pow(b2x,2)*b3x*c2x*c2y*c3y + 2*pow(b2y,2)*b3x*c2x*c2y*c3y + 4*b2x*b2y*b3y*c2x*c2y*c3y - 
		    2*b2x*b2y*b3x*pow(c2y,2)*c3y - pow(b2x,2)*b3y*pow(c2y,2)*c3y + pow(b2y,2)*b3y*pow(c2y,2)*c3y);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,b2w*b2w*b3w*a2w*a2w*a3w);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,c2w*c2w*c3w*a2w*a2w*a3w);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,c2w*c2w*c3w*b2w*b2w*b3w);

  qanal[anal].qBA[bin][k]->Fill(0.,pow(a2x,2)*a3x*pow(b2x,2)*b3x - pow(a2y,2)*a3x*pow(b2x,2)*b3x - 2*a2x*a2y*a3y*pow(b2x,2)*b3x + 4*a2x*a2y*a3x*b2x*b2y*b3x + 2*pow(a2x,2)*a3y*b2x*b2y*b3x - 
		    2*pow(a2y,2)*a3y*b2x*b2y*b3x - pow(a2x,2)*a3x*pow(b2y,2)*b3x + pow(a2y,2)*a3x*pow(b2y,2)*b3x + 2*a2x*a2y*a3y*pow(b2y,2)*b3x + 2*a2x*a2y*a3x*pow(b2x,2)*b3y + 
		    pow(a2x,2)*a3y*pow(b2x,2)*b3y - pow(a2y,2)*a3y*pow(b2x,2)*b3y - 2*pow(a2x,2)*a3x*b2x*b2y*b3y + 2*pow(a2y,2)*a3x*b2x*b2y*b3y + 4*a2x*a2y*a3y*b2x*b2y*b3y - 
		    2*a2x*a2y*a3x*pow(b2y,2)*b3y - pow(a2x,2)*a3y*pow(b2y,2)*b3y + pow(a2y,2)*a3y*pow(b2y,2)*b3y);
  qanal[anal].qCA[bin][k]->Fill(0.,pow(a2x,2)*a3x*pow(c2x,2)*c3x - pow(a2y,2)*a3x*pow(c2x,2)*c3x - 2*a2x*a2y*a3y*pow(c2x,2)*c3x + 4*a2x*a2y*a3x*c2x*c2y*c3x + 2*pow(a2x,2)*a3y*c2x*c2y*c3x - 
		    2*pow(a2y,2)*a3y*c2x*c2y*c3x - pow(a2x,2)*a3x*pow(c2y,2)*c3x + pow(a2y,2)*a3x*pow(c2y,2)*c3x + 2*a2x*a2y*a3y*pow(c2y,2)*c3x + 2*a2x*a2y*a3x*pow(c2x,2)*c3y + 
		    pow(a2x,2)*a3y*pow(c2x,2)*c3y - pow(a2y,2)*a3y*pow(c2x,2)*c3y - 2*pow(a2x,2)*a3x*c2x*c2y*c3y + 2*pow(a2y,2)*a3x*c2x*c2y*c3y + 4*a2x*a2y*a3y*c2x*c2y*c3y - 
		    2*a2x*a2y*a3x*pow(c2y,2)*c3y - pow(a2x,2)*a3y*pow(c2y,2)*c3y + pow(a2y,2)*a3y*pow(c2y,2)*c3y);
  qanal[anal].qCB[bin][k]->Fill(0.,pow(b2x,2)*b3x*pow(c2x,2)*c3x - pow(b2y,2)*b3x*pow(c2x,2)*c3x - 2*b2x*b2y*b3y*pow(c2x,2)*c3x + 4*b2x*b2y*b3x*c2x*c2y*c3x + 2*pow(b2x,2)*b3y*c2x*c2y*c3x - 
		    2*pow(b2y,2)*b3y*c2x*c2y*c3x - pow(b2x,2)*b3x*pow(c2y,2)*c3x + pow(b2y,2)*b3x*pow(c2y,2)*c3x + 2*b2x*b2y*b3y*pow(c2y,2)*c3x + 2*b2x*b2y*b3x*pow(c2x,2)*c3y + 
		    pow(b2x,2)*b3y*pow(c2x,2)*c3y - pow(b2y,2)*b3y*pow(c2x,2)*c3y - 2*pow(b2x,2)*b3x*c2x*c2y*c3y + 2*pow(b2y,2)*b3x*c2x*c2y*c3y + 4*b2x*b2y*b3y*c2x*c2y*c3y - 
		    2*b2x*b2y*b3x*pow(c2y,2)*c3y - pow(b2x,2)*b3y*pow(c2y,2)*c3y + pow(b2y,2)*b3y*pow(c2y,2)*c3y);
  qanal[anal].qBAcnt[bin][k]->Fill(0.,b2w*b2w*b3w*a2w*a2w*a3w);
  qanal[anal].qCAcnt[bin][k]->Fill(0.,c2w*c2w*c3w*a2w*a2w*a3w);
  qanal[anal].qCBcnt[bin][k]->Fill(0.,c2w*c2w*c3w*b2w*b2w*b3w);


}
