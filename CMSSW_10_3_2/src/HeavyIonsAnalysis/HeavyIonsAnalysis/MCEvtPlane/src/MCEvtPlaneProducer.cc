// -*- C++ -*-
//
// Package:    MCEvtPlaneProducer
// Class:      MCEvtPlaneProducer
// 
/**\class MCEvtPlaneProducer MCEvtPlaneProducer.cc RecoHI/MCEvtPlaneProducer/src/MCEvtPlaneProducer.cc
   
Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Sergey Petrushanko
//         Created:  Fri Jul 11 10:05:00 2008
//
//

// system include files
#include <memory>
#include <iostream>
#include <time.h>
#include <cmath>
#include "TMath.h"
#include "TF1.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Candidate/interface/Candidate.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"


#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CastorReco/interface/CastorTower.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include <cstdlib>
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

using namespace std;
using namespace hi;
using namespace edm;
using namespace reco;
//
// class decleration
//

namespace hi {
  class GenPlane {
  public:
    GenPlane(string name,double etaminval1,double etamaxval1,double etaminval2,double etamaxval2,int orderval){
      epname=name;
      etamin1=etaminval1;
      etamax1=etamaxval1;
      etamin2=etaminval2;
      etamax2=etamaxval2;
      sumsin=0;
      sumcos=0;
      sumsinNoWgt=0;
      sumcosNoWgt=0;
      
      mult = 0;
      order = (double) orderval;
    }
    ~GenPlane(){;}
    void addParticle(double w, double PtOrEt, double s, double c, double eta) {
      if((eta>=etamin1 && eta<etamax1) || 
	 (etamin2!= etamax2 && eta>=etamin2 && eta<etamax2 )) {
	sumsin+=w*s;
	sumcos+=w*c;
	sumsinNoWgt+=s;
	sumcosNoWgt+=c;
	
	sumw+=fabs(w);
	sumw2+=w*w;
	sumPtOrEt+=PtOrEt;
	sumPtOrEt2+=PtOrEt*PtOrEt;
	++mult;
      }
    }
    
    double getAngle(double &ang, double &sv, double &cv, double &svNoWgt, double &cvNoWgt,  double &w, double &w2, double &PtOrEt, double &PtOrEt2, uint &epmult){
      ang = -10;
      sv = 0;
      cv = 0;
      sv = sumsin;
      cv = sumcos;
      svNoWgt = sumsinNoWgt;
      cvNoWgt = sumcosNoWgt;
      w = sumw;
      w2 = sumw2;
      PtOrEt = sumPtOrEt;
      PtOrEt2 = sumPtOrEt2;
      epmult = mult;
      double q = sv*sv+cv*cv;
      if(q>0) ang = atan2(sv,cv)/order;
      return ang;
    }
    void reset() {
      sumsin=0;
      sumcos=0;
      sumsinNoWgt = 0;
      sumcosNoWgt = 0;
      sumw = 0;
      sumw2 = 0;
      mult = 0;
      sumPtOrEt = 0;
      sumPtOrEt2 = 0;
    }
  private:
    string epname;
    double etamin1;
    double etamax1;
    
    double etamin2;
    double etamax2;
    double sumsin;
    double sumcos;
    double sumsinNoWgt;
    double sumcosNoWgt;
    uint mult;
    double sumw;
    double sumw2;
    double sumPtOrEt;
    double sumPtOrEt2;
    double order;
  };
}  

class MCEvtPlaneProducer : public edm::stream::EDProducer<> {
public:
  explicit MCEvtPlaneProducer(const edm::ParameterSet&);
  int bToCent(int b);
  ~MCEvtPlaneProducer();

private:
  GenPlane *rp[NumEPNames];

  void produce(edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------

  edm::InputTag					trackEta_;
  edm::InputTag					trackPhi_;
  edm::InputTag					trackPt_;
  edm::InputTag					trackWeight_;
  edm::InputTag					vertexZ_;
  edm::InputTag					chargeTag_;
  edm::InputTag RP_;
  edm::InputTag Ncoll_;
  edm::InputTag Npart_;
  bool loadDB_;
  edm::InputTag b_;
  std::string centralityVariable_;
  std::string centralityLabel_;
  std::string centralityMC_;
  int FlatOrder_;
  int NumFlatBins_;
  int flatnvtxbins_;
  double flatminvtx_;
  double flatdelvtx_;
  vector<double> bToCentCalib_;
  HiEvtPlaneFlatten * flat[NumEPNames];
  // int ntrkval;
  // int NtrkToBin(int ntrk){
  //   for(int i = 0; i<=ntrkbins; i++) {
  //     if(ntrk < trkBins[i]) return i;
  //   }
  //   return ntrkbins;
  // }
};

MCEvtPlaneProducer::MCEvtPlaneProducer(const edm::ParameterSet& iConfig):
  trackEta_( iConfig.getUntrackedParameter<edm::InputTag>("trackEta") ),
  trackPhi_( iConfig.getUntrackedParameter<edm::InputTag>("trackPhi") ),
  trackPt_( iConfig.getUntrackedParameter<edm::InputTag>("trackPt") ),
  trackWeight_( iConfig.getUntrackedParameter<edm::InputTag>("trackWeight") ),
  vertexZ_( iConfig.getUntrackedParameter<edm::InputTag>("vertexZ") ),
  chargeTag_( iConfig.getUntrackedParameter<edm::InputTag>("trackCharge") ),
  RP_( iConfig.getUntrackedParameter<edm::InputTag>("eventRP") ),
  Ncoll_( iConfig.getUntrackedParameter<edm::InputTag>("eventNcoll") ),
  Npart_( iConfig.getUntrackedParameter<edm::InputTag>("eventNpart") ),
  loadDB_ ( iConfig.getParameter<bool>("loadDB") ),
  b_( iConfig.getUntrackedParameter<edm::InputTag>("eventb") ),
  FlatOrder_ ( iConfig.getUntrackedParameter<int>("FlatOrder") ),
  NumFlatBins_ ( iConfig.getUntrackedParameter<int>("NumFlatBins") ),
  flatnvtxbins_ ( iConfig.getParameter<int>("flatnvtxbins") ),
  flatminvtx_ ( iConfig.getParameter<double>("flatminvtx") ),
  flatdelvtx_ ( iConfig.getParameter<double>("flatdelvtx") ),
  bToCentCalib_ ( iConfig.getParameter<std::vector<double>>("bToCentCalib") )
{
  consumes<double>(RP_);
  consumes<double>(Ncoll_);
  consumes<double>(Npart_);
  consumes<double>(b_);
  consumes<std::vector<double> >(trackEta_);
  consumes<std::vector<double> >(trackPhi_);
  consumes<std::vector<double> >(trackPt_);
  consumes<std::vector<double> >(trackWeight_);
  consumes<std::vector<double> >(vertexZ_);
  consumes<std::vector<double> >(chargeTag_);
  
  //nCentBins_ = 200.;


  produces<reco::EvtPlaneCollection>();
  produces<reco::VertexCollection>();
  produces<reco::TrackCollection>();
  produces<reco::Centrality>();
  produces<int>("centralityBin");
  produces<double>("b");
  produces<double>("rp");
  for(int i = 0; i<NumEPNames; i++ ) {
    rp[i] = new GenPlane(EPNames[i].data(),EPEtaMin1[i],EPEtaMax1[i],EPEtaMin2[i],EPEtaMax2[i],EPOrder[i]);
  }
  for(int i = 0; i<NumEPNames; i++) {
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->init(FlatOrder_,NumFlatBins_,flatnvtxbins_,flatminvtx_,flatdelvtx_,EPNames[i],EPOrder[i]);
  }
}


MCEvtPlaneProducer::~MCEvtPlaneProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MCEvtPlaneProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  //
  //Get flattening parameter file.  
  //
  if ( loadDB_ ) {
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    LoadEPDB db(flatparmsDB_,flat);
    if(!db.IsSuccess()) {
      loadDB_ = kFALSE;
    }
  }
  
  Handle<std::vector<double> >	fEta;
  Handle<std::vector<double> >	fPhi;
  Handle<std::vector<double> >	fPt;
  Handle<std::vector<double> >	fWeight;
  Handle<std::vector<double> >	fVz;
  Handle<std::vector<double> >  fCharge;
  
  iEvent.getByLabel(trackEta_,	fEta);
  iEvent.getByLabel(trackPhi_,	fPhi);
  iEvent.getByLabel(trackPt_,	fPt);
  iEvent.getByLabel(trackWeight_, fWeight);
  iEvent.getByLabel(vertexZ_, 	fVz);
  iEvent.getByLabel(chargeTag_, 	fCharge);
  
  Handle<double> fRP;
  Handle<double> fNcoll;
  Handle<double> fNpart;
  Handle<double> fb;
  iEvent.getByLabel(RP_, fRP);
  iEvent.getByLabel(Ncoll_, fNcoll);
  iEvent.getByLabel(Npart_, fNpart);
  iEvent.getByLabel(b_, fb);
  for(int i = 0; i<NumEPNames; i++) rp[i]->reset();
  int bin = bToCent(500.* (*fb));
  if(bin<0) bin = 0;
  if(bin>=200) bin=199;

  iEvent.put(std::make_unique<int>(bin),"centralityBin");
  iEvent.put(std::make_unique<double>(*fb),"b");
  iEvent.put(std::make_unique<double>(*fRP-TMath::Pi()),"rp");
  //Create Vertex
  const double chi2 = 20.0;
  const double ndof = 10;
  const size_t size = 0; // can't test refs at the moment
  const double x = 1.0, y = 2.0, z = 0.0;
  const float ee[ reco::Vertex::Error::kSize ] = { 1.1, 1.2, 1.3, 
						   2.2, 2.3, 
						   3.3 };
  reco::Vertex::Error err;
  int k = 0;
  for( int i = 0; i < reco::Vertex::Error::kRows; ++i )
    for( int j = i; j < reco::Vertex::Error::kCols; ++j )
      err( i, j ) = ee[ k++ ];
  
  reco::Vertex v( reco::Vertex::Point( x, y, z ), err, chi2, ndof, size );
  auto vtxOutput = std::make_unique<VertexCollection>();
  vtxOutput->push_back(v);
  iEvent.put(std::move(vtxOutput));
  
  //Tracking part
  double track_eta;
  double track_phi;
  double track_pt;
  double track_charge;
  //double track_vz;
  double minpt_ = -1;
  double maxpt_ = -1;
  const reco::Vertex::Point vt(0,0,0);
  float vzr_sell = 0;
  typedef math::Error<5>::type CovarianceMatrix;
  const CovarianceMatrix cov;

  auto trkOutput = std::make_unique<TrackCollection>();
  int icount = 0;
  for(int j = 0; j<int(fEta->size()); j++) {
    track_eta = (*fEta)[j];
    track_phi = (*fPhi)[j];
    track_pt = (*fPt)[j];
    track_charge = (*fCharge)[j];
    ++icount;
    //std::cout<<"eta,phi,pt,charge: "<<track_eta<<"\t"<<track_phi<<"\t"<<track_pt<<"\t"<<track_charge<<std::endl;
    double theta = 2*TMath::ATan(TMath::Exp(-track_eta));

    const TrackExtra::Vector  momentum(track_pt*cos(track_phi),track_pt*sin(track_phi),track_pt/TMath::Tan(theta));
    //momentum[1]=track_pt*sin(track_phi);
    //track_vz = (*fVz)[j];
    reco::Track t(0.,0.,vt,momentum,(int)track_charge,cov,reco::TrackBase::initialStep, reco::TrackBase::highPurity);
    trkOutput->push_back(t);
    double minpt = minpt_;
    double maxpt = maxpt_;
    for(int i = 0; i<NumEPNames; i++) {
      if(minpt_<0) minpt = minTransverse[i];
      if(maxpt_<0) maxpt = maxTransverse[i];
      if(track_pt<minpt) continue;
      if(track_pt>maxpt) continue;
      //      if(EPDet[i]==Tracker) {
	double w = track_pt;
	if(w>2.5) w=2.0;   //v2 starts decreasing above ~2.5 GeV/c
	if(EPOrder[i]==1) {
	 if(MomConsWeight[i][0]=='y' && loadDB_) {
	  w = flat[i]->getW(track_pt, vzr_sell, bin);
	  }
	 //if(track_eta<0 && MomConsWeight[i][0]=='n') w=-w;
	}
	rp[i]->addParticle(w,track_pt,sin(EPOrder[i]*track_phi),cos(EPOrder[i]*track_phi),track_eta);
	//      }
    }
  }
  auto evtplaneOutput = std::make_unique<EvtPlaneCollection>();
  double ang=-10;
  double sv = 0;
  double cv = 0;
  double svNoWgt = 0;
  double cvNoWgt = 0;
  
  double wv = 0;
  double wv2 = 0;
  double pe = 0;
  double pe2 = 0;
  uint epmult = 0;
  
  for(int i = 0; i<NumEPNames; i++) {
    rp[i]->getAngle(ang,sv,cv,svNoWgt, cvNoWgt, wv,wv2,pe,pe2,epmult);
    EvtPlane *ep = new EvtPlane(i,2,ang,sv,cv,wv,wv2,pe,pe2,epmult) ;
    ep->addLevel(0,ang,sv,cv);
    ep->addLevel(1,ang,sv,cv);
    evtplaneOutput->push_back( *ep );
    evtplaneOutput->back().addLevel(3, 0., svNoWgt, cvNoWgt);

  }
  
  iEvent.put(std::move(evtplaneOutput));
  iEvent.put(std::move(trkOutput));
  auto centOutput = std::make_unique<Centrality>();
  iEvent.put(std::move(centOutput));
}
int MCEvtPlaneProducer::bToCent(int b){
  for(int i=199; i>0; i--) {
    if(b<bToCentCalib_.at(i)) continue;
    return i;
  }
  return 0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MCEvtPlaneProducer);
