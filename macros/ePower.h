#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "Math/DistFunc.h"
#include "Math/SpecFuncMathMore.h"
#include "TF1.h"

#include <iostream>

///////////////////////////////////////////////////////////////////////////

double ePower(double vn, double alpha, double kn, double e0){
  //cout<<vn<<"\t"<<alpha<<"\t"<<kn<<"\t"<<e0<<endl;
  //if(alpha<=0||kn<=0||e0<=0||isnan(alpha)||isnan(kn)||isnan(e0)) return 1e5;
  double hyper = ROOT::Math::hyperg(0.5,1 + 2*alpha,1,(2*e0*vn)/(kn + e0*vn));
  if(isnan(hyper)||hyper==0) cout<<1+2*alpha<<"\t"<<(2*e0*vn)/(kn + e0*vn)<<endl;
  double res = (2*alpha*pow(1 - pow(e0,2),0.5 + alpha)*vn*pow(1 + (e0*vn)/kn,-1 - 2*alpha)*
		pow(1 - pow(vn,2)/pow(kn,2),-1 + alpha)*hyper)/pow(kn,2);
  //cout<<"hyperg: "<<res<<"\t"<<hyper<<"\t"<<vn<<"\t"<<alpha<<"\t"<<kn<<"\t"<<e0<<"\t"<<(2*e0*vn)/(kn + e0*vn)<<endl;
  
  return res;
}

Double_t ePowerWrap(Double_t *x, Double_t *par) {
  double xx = x[0];
  double alpha = par[0];
  double kn = par[1];
  double e0 = par[2];
  double scale = par[3];
  double f = scale*ePower(xx,alpha,kn,e0);
  if(f<1e-6 || isnan(f)) return 1e-6;
  return f;
}

TF1 * ePowerFunc = new TF1("ePowerWrap",ePowerWrap,0,1,4);

///////////////////////////////////////////////////////////////////////////

double stuT(double t, double nu, double mu, double delta){
  double num = ROOT::Math::tgamma((nu+1)/2.);
  double denom = ROOT::Math::tgamma(nu/2.) * delta * TMath::Sqrt(nu*TMath::Pi());
  double w = nu*pow(delta,2);
  double res=(num/denom)*pow(1 + pow(t-mu,2.)/w,-(nu+1)/2);
  return res;

}

double stuTWrap(double * x, double * par){
  double xx = x[0];
  double nu = par[0];
  double mu = par[1];
  double delta = par[2];
  double scale = par[3];
  return scale*stuT(xx,nu,mu,delta);
}

double stuTcore(double phi, double vn, double vobs, double delta, double nu){
  double num = vobs*ROOT::Math::tgamma((nu+1)/2.);
  double denom = ROOT::Math::tgamma(nu/2.) * delta * TMath::Sqrt(nu*TMath::Pi());
  double w = nu*pow(delta,2);
  double arg = pow(1 + (pow(vobs,2)+pow(vn,2)-2*vobs*vn*TMath::Cos(phi))/w,-(nu+1)/2);
  double res = num*arg/denom;
  return res;
}


double stuT2d(double vn, double vobs, double delta, double nu){
  double dphi = 2*TMath::Pi()/5000.;
  double sum = 0;
  double phi = 0;
  while (phi<2*TMath::Pi()) {
    double newterm = stuTcore(phi,vn,vobs,delta,nu);
    sum+=newterm;
    phi+=dphi;
  }
  return sum*dphi;
}

double stuT2dWrap(double * x, double * par){
  double xx = x[0];
  double vobs = par[0];
  double delta = par[1];
  double nu = par[2];
  double scale = par[3];
  return scale*stuT2d(xx,vobs,delta,nu);
}

TF1 * stuT2dFunc = new TF1("stuT2dFunc",stuT2dWrap,0,1,4);

///////////////////////////////////////////////////////////////////////////

double gaus2d(double vn, double vobs, double delta) {
  double res = (vobs*ROOT::Math::cyl_bessel_i(0,(vn*vobs)/pow(delta,2)))/(pow(delta,2)*pow((double)TMath::E(),(pow(vn,2) + pow(vobs,2))/(2.*pow(delta,2))));
  return res;
}

double gaus2dWrap(double * x, double * par) {
  double vobs = x[0];
  double vn = par[0];
  double delta = par[1];
  double scale = par[2];
  return scale*gaus2d(vn,vobs,delta);
}

TF1 * gaus2dFunc = new TF1("gaus2dFunc",gaus2dWrap,0,1,3);

///////////////////////////////////////////////////////////////////////////

double vnGaus(double vobs, double delta=0.15, double alpha=23.2, double kn=0.3380, double e0=0.2641){
  double dvn = 0.001;
  double vnsum = 0;
  double vn = dvn;
  double v = 10;
  while(vn<1 && v>1e-10){
    double g = gaus2d(vn,vobs,delta);
    v = ePower(vn,alpha,kn,e0)*g;
    vnsum+=v;
    vn+=dvn;
  }
  //  cout<<"vnsum: "<<vnsum<<endl;
  return vnsum*dvn;
}

double vnGausWrap(double * x, double * par){
  double vobs = x[0];
  double delta = par[0];
  double alpha = par[1];
  double kn = par[2];
  double e0 = par[3];
  double scale = par[4];
  if(delta==0) return scale*ePower(vobs,alpha,kn,e0);
  return scale*vnGaus(vobs,delta,alpha,kn,e0);
}

TF1 * vnGausFunc = new TF1("vnGausFunc",vnGausWrap,0,1,5);

///////////////////////////////////////////////////////////////////////////

double vnStuT(double vobs, double delta=0.02, double nu=10, double alpha=23.2, double kn=0.338, double e0=0.2641){
  double dvn = 0.001;
  double vn = 0.0;
  double vnsum = 0;
  while(vn<1){
    double v = ePower(vn,alpha,kn,e0)*stuT2d(vn,vobs,delta,nu);
    vnsum+=v;
    vn+=dvn;
  }
  return vnsum*dvn;
}

double vnStuTWrap(double *x, double *par){
  double vobs = x[0];
  double delta = par[0];
  double nu = par[1];
  double alpha = par[2];
  double kn = par[3];
  double e0 = par[4];
  return vnStuT(vobs,delta,nu,alpha,kn,e0);
}

TF1 * vnStuTFunc = new TF1("vnStuTFunc",vnStuTWrap,0,1,5);
