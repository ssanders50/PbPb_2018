// -*- C++ -*-
//
// Package:    VNAnalyzer
// Class:      VNAnalyzer
// 

// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HeavyIonsAnalysis/VNAnalysis/interface/TrackEfficiency.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <time.h>
#include <cstdlib>
	
#include <vector>
using std::vector;
using std::rand;
using namespace std;
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/TrackQuality.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"
using namespace hi;
using namespace edm;


static const int ntrkbins = 15;
static const  float trkBins[]={0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100};

static const int netabinsDefault = 12;
static const float etabinsDefault[]={-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

static const int ncentbins = 13;
static const  float centbins[]={0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100};

static const int nptbins = 28;
static const float ptbins[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0,
			      2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  7.0, 8.0, 10., 12.0, 14.0, 16.0,  20.0, 26.0, 35.0, 45.0, 60.0, 80.0, 100., 200.};

//
// class declaration
//

class VNAnalyzer : public edm::EDAnalyzer {
public:
  explicit VNAnalyzer(const edm::ParameterSet&);
  ~VNAnalyzer();
      
private:
  int NtrkToBin(int ntrk){
    for(int i = 0; i<=ntrkbins; i++) {
      if(ntrk < trkBins[i]) return i;
    }
    return ntrkbins;
  }
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool CaloMatch(const reco::Track & track, const edm::Event & iEvent, unsigned int idx);
  // ----------member data ---------------------------
  int eporder_;
  int mx;
  std::string centralityVariable_;
  std::string centralityLabel_;
  std::string centralityMC_;
  edm::InputTag centralityBinTag_;
  edm::EDGetTokenT<int> centralityBinToken;
  edm::Handle<int> cbin_;
  edm::EDGetTokenT<int> tag_;

  edm::InputTag bTag_;
  edm::EDGetTokenT<double> bToken;
  edm::Handle<double> b_;
  double b;

  edm::InputTag RPTag_;
  edm::EDGetTokenT<double> RPToken;
  edm::Handle<double> RP_;
  double rp;

  edm::InputTag centralityTag_;
  edm::EDGetTokenT<reco::Centrality> centralityToken;
  edm::Handle<reco::Centrality> centrality_;

  edm::InputTag pfTag_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfToken_;

  edm::InputTag vertexTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;

  edm::InputTag inputPlanesTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanes_;

  edm::Service<TFileService> fs;
  TFile * frecenter=0;
  string offsetFileName;
 
  double caloCentRef_;
  double caloCentRefWidth_;
  int caloCentRefMinBin_;
  int caloCentRefMaxBin_;

  bool useNtrk_;

  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH1D * hcent;
  TH1D * hb;
  TH1D * hRP;
  TH1D * hvtx;
  TH1D * hvtxRaw;
  TH1D * hcentbins;
  TH1D * hcentres;
  TH1D * hptNtrk;
  TH1D * hptNtrkGood;
  TH1I * hNtrkRet;
  //  TH2D * hEff[ntrkbins];
  double centval;
  int ntrkval;
  int Noff;
  double vtx;

  double reso_;
  bool bCaloMatching_;
  int nvtx_;
  double minvz_;
  double maxvz_;
  double dzdzerror_pix_;
  double chi2_;

  double dzdzerror_;
  double d0d0error_;
  double pterror_;
  TrackQuality * trackq;

  Double_t epang[NumEPNames];
  Double_t eporig[NumEPNames];
  Double_t epsin[NumEPNames];
  Double_t epcos[NumEPNames];

  Double_t qx[NumEPNames];
  Double_t qy[NumEPNames];
  Double_t q[NumEPNames];
  Double_t epmult[NumEPNames];
  Double_t sumw[NumEPNames];
  Double_t sumw2[NumEPNames];
  Double_t vn[NumEPNames];

  Double_t rescor[NumEPNames];
  Double_t rescorErr[NumEPNames];
  TH1D * hPsi[NumEPNames];
  TH1D * hPsiOffset[NumEPNames];
  TH1D * hPsiFlat[NumEPNames];


  unsigned int runno_;

  TH1D * hNtrk;
  TH1D * hNoff;
  int nEtaBins;
  TH1I * hrun;
  string rpnames[NumEPNames];
  string effTable_;
  TTree * tree;
  TrackEfficiency *teff;
  int FlatOrder_;
  int NumFlatBins_;
  int CentBinCompression_;
  int flatnvtxbins_;
  double flatminvtx_;
  double flatdelvtx_;
  bool genMC_;
  int Noffmin_;
  int Noffmax_;
  TH2F * qxtrk[7];
  TH2F * qytrk[7];
  TH2F * qcnt;
  TH2F * avpt;
  TH2D * res[7][ntrkbins];
  TH2D * resw[7][ntrkbins];
  TH2D * resep[7][ntrkbins];
  TH2D * rescnt[7][ntrkbins];
  TH2D * res12[ntrkbins];
  TH2D * resw12[ntrkbins];
  TH2D * resep12[ntrkbins];
  TH2D * rescnt12[ntrkbins];
  TH2D * ptspec[ntrkbins];  
  HiEvtPlaneFlatten * flat[NumEPNames];
  bool loadDB_;
  bool useNtrkBins_; 
  bool bypassCentrality_;
  bool FirstEvent_;
  bool makeTree_;
  bool Recenter_;
  int minrun_;
  int maxrun_;
  struct offsets {
    TH2D * wqxtrkRef[7][40];
    TH2D * wqytrkRef[7][40];
  } offsets[12];
  int nCentBins_ = 1;
  int ntrack;


  TH2D * hTemplate;
  //===================================

  
  TH2D * Eff_0_5;
  TH2D * Eff_5_10;
  TH2D * Eff_10_30;
  TH2D * Eff_30_50;
  TH2D * Eff_50_100;

  TRandom * ran;
  
  //===================================
  
  int getNoff(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    int N = 0;
    using namespace edm;
    using namespace reco;
    iEvent.getByToken(vertexToken,vertex_);
    VertexCollection recoVertices = *vertex_;
    int primaryvtx = 0;
    TrackQuality * tracknoff = new TrackQuality(trackTag_, 3., 3., 0.1, 10., 40.);
    if(tracknoff->is_ppReco()) {
      if ( recoVertices.size() > 100 ) return -1;
      sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
	  if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
	  return a.tracksSize() > b.tracksSize();
	});
    } else {
      sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
	  return a.tracksSize() > b.tracksSize();
	});    
    }
    
    vtx = recoVertices[primaryvtx].z();
    if (fabs(vtx) < minvz_ || fabs(vtx) > maxvz_) {
      return -1;
    }        
    
    iEvent.getByLabel(trackTag_,trackCollection_);
    
    for(TrackCollection::const_iterator itTrack = trackCollection_->begin(); itTrack != trackCollection_->end(); ++itTrack) { 
      if(!tracknoff->isGood(itTrack,recoVertices)) continue;
      
      if( itTrack->pt() < 0.4 ) continue;
      ++N;
    }
    return N;
  }
  
  
  int fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup, int bin)
  {
    int Ntrk = 0;
    using namespace edm;
    using namespace reco;
    for(int i = 0; i<7; i++) {
      qxtrk[i]->Reset();
      qytrk[i]->Reset();
    }
    qcnt->Reset();
    avpt->Reset();
    iEvent.getByToken(vertexToken,vertex_);
    VertexCollection recoVertices = *vertex_;
    int primaryvtx = 0;
    if(trackq->is_ppReco()) {
      if ( recoVertices.size() > 100 ) return -1;
      sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
	  if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
	  return a.tracksSize() > b.tracksSize();
	});
    } else {
      sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
	  return a.tracksSize() > b.tracksSize();
	});    
    }
          
    vtx = recoVertices[primaryvtx].z();
    hvtxRaw->Fill(vtx);

    if (fabs(vtx) < minvz_ || fabs(vtx) > maxvz_) {
      return -1;
    }        
    
    iEvent.getByLabel(trackTag_,trackCollection_);
    int k = (vtx-flatminvtx_)/flatdelvtx_;
    for(TrackCollection::const_iterator itTrack = trackCollection_->begin(); itTrack != trackCollection_->end(); ++itTrack) { 
      if(!trackq->isGood(itTrack,recoVertices)) continue;
      double eta = itTrack->eta();
      double pt = itTrack->pt();
      double phi = itTrack->phi();
      int ipt = qxtrk[0]->GetXaxis()->FindBin(itTrack->pt());
      int ieta = qxtrk[0]->GetYaxis()->FindBin(itTrack->eta());
      double eff = 1.;
      // if(effTable_!="NULL") {
      // 	int ieffpt = Eff_0_5->GetYaxis()->FindBin(itTrack->pt());
      // 	int ieffeta = Eff_0_5->GetXaxis()->FindBin(itTrack->eta());
      // 	if(centval<5) eff = Eff_0_5->GetBinContent(ieffeta,ieffpt);
      // 	else if (centval < 10) eff = Eff_5_10->GetBinContent(ieffeta,ieffpt);
      // 	else if (centval < 30) eff = Eff_10_30->GetBinContent(ieffeta,ieffpt);
      // 	else if (centval < 50) eff = Eff_30_50->GetBinContent(ieffeta,ieffpt);
      // 	else eff = Eff_50_100->GetBinContent(ieffeta,ieffpt);	
      // 	if(eff == 0) eff = 1;
      // 	eff=1/eff;
      // }
      if(k>=0&&k<flatnvtxbins_) {
	for(int iorder = 1; iorder <=7; iorder++) {
	  qxtrk[iorder-1]->Fill(pt, eta, eff*(TMath::Cos(iorder*phi) - offsets[k].wqxtrkRef[iorder-1][bin]->GetBinContent(ipt,ieta)));
	  qytrk[iorder-1]->Fill(pt, eta, eff*(TMath::Sin(iorder*phi) - offsets[k].wqytrkRef[iorder-1][bin]->GetBinContent(ipt,ieta)));
	}
	qcnt->Fill(itTrack->pt(), itTrack->eta(), eff);
	ptspec[bin]->Fill(pt,eta);
	avpt->Fill(itTrack->pt(), itTrack->eta(), eff*itTrack->pt());
      }
      if( itTrack->pt() < 0.2 ) continue;
      ++Ntrk;
    }
    return Ntrk;
  }
};


//
// constructors and destructor
//
VNAnalyzer::VNAnalyzer(const edm::ParameterSet& iConfig):runno_(0)
							 
{
  ran = new TRandom();
  ran->SetSeed(0);
  runno_ = 0;
  loadDB_ = kTRUE;
  FirstEvent_ = kTRUE;
  for(int i = 0; i<NumEPNames; i++) {
    epang[i] = -10;
    eporig[i] = -10;
    epsin[i] = 0;
    epcos[i] = 0;
    qx[i] = 0;
    qy[i] = 0;
    q[i] = 0;
    epmult[i] = 0;
    rescor[i] = 0;
    rescorErr[i] = 0;
  }
  
  centralityVariable_ = iConfig.getParameter<std::string>("centralityVariable");
  if(iConfig.exists("nonDefaultGlauberModel")){
    centralityMC_ = iConfig.getParameter<std::string>("nonDefaultGlauberModel");
  }
  centralityLabel_ = centralityVariable_+centralityMC_;
  
  centralityBinTag_ = iConfig.getParameter<edm::InputTag>("centralityBinTag_");
  centralityBinToken = consumes<int>(centralityBinTag_);
  
  bTag_ = iConfig.getParameter<edm::InputTag>("bTag_");
  bToken = consumes<double>(bTag_);
  
  RPTag_ = iConfig.getParameter<edm::InputTag>("RPTag_");
  RPToken = consumes<double>(RPTag_);
  
  centralityTag_ = iConfig.getParameter<edm::InputTag>("centralityTag_");
  centralityToken = consumes<reco::Centrality>(centralityTag_);
  if(centralityToken.isUninitialized()) {
    std::cout<<"centralityToken is uninitialized."<<std::endl;
  }
  vertexTag_  = iConfig.getParameter<edm::InputTag>("vertexTag_");
  vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);
  if(vertexToken.isUninitialized()) {
    std::cout<<"vertexToken is uninitialized."<<std::endl;
  }
  
  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
  trackToken = consumes<reco::TrackCollection>(trackTag_);
  if(trackToken.isUninitialized()) {
    std::cout<<"trackToken is uninitialized."<<std::endl;
  }
  useNtrk_ = iConfig.getUntrackedParameter<bool>("useNtrk",false);
  if(useNtrk_) {
    NumFlatBins_ = ntrkbins;
    CentBinCompression_ = 1;
  }
  
  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesTag_");
  inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);
  if(inputPlanesToken.isUninitialized()) {
    std::cout<<"inputPlanesToken is uninitialized."<<std::endl;
  }
  tag_ = consumes<int>(iConfig.getParameter<edm::InputTag>("BinLabel"));
  
  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatBins_ = iConfig.getUntrackedParameter<int>("NumFlatBins_",20);
  caloCentRef_ = iConfig.getUntrackedParameter<double>("caloCentRef_",80.);
  caloCentRefWidth_ = iConfig.getUntrackedParameter<double>("caloCentRefWidth_",5.);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
  flatnvtxbins_ = iConfig.getParameter<int>("flatnvtxbins") ;
  flatminvtx_ = iConfig.getParameter<double>("flatminvtx") ;
  flatdelvtx_ = iConfig.getParameter<double>("flatdelvtx") ;
  genMC_ = iConfig.getUntrackedParameter<bool>("genMC",false);
  
  Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
  Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 50000);	
  minrun_ = iConfig.getUntrackedParameter<int>("minrun_", 0);
  maxrun_ = iConfig.getUntrackedParameter<int>("maxrun_", 50000);	
  effTable_ = iConfig.getParameter<std::string>("effTable_");
  bCaloMatching_ = iConfig.getUntrackedParameter<bool>("bCaloMaching", false);
  Recenter_ = iConfig.getUntrackedParameter<bool>("Recenter",true);
  makeTree_ = iConfig.getUntrackedParameter<bool>("makeTree_",false);
  
  nvtx_ = iConfig.getUntrackedParameter<int>("nvtx_", 100);
  reso_ = iConfig.getUntrackedParameter<double>("reso", 0.2);
  if(reso_<0.01) bCaloMatching_ = false;
  if(bCaloMatching_) {
    pfTag_ = iConfig.getUntrackedParameter<edm::InputTag>("pfTag");
    pfToken_ = consumes<reco::PFCandidateCollection>(pfTag_);
  }
  dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror_", 3.);
  d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error_", 3.);
  pterror_ = iConfig.getUntrackedParameter<double>("pterror_",0.1);
  dzdzerror_pix_ = iConfig.getUntrackedParameter<double>("dzdzerror_pix_") ;
  chi2_  = iConfig.getUntrackedParameter<double>("chi2_") ;
  
  trackq = new TrackQuality(trackTag_, dzdzerror_, d0d0error_, pterror_, dzdzerror_pix_, chi2_);  
  teff = 0;
  if(effTable_!="NULL") teff = new TrackEfficiency(effTable_.data());
  minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
  maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
  offsetFileName = iConfig.getUntrackedParameter<std::string>("offsetFile");
  if(!genMC_ && offsetFileName.length()!=std::string::npos) {
    frecenter = new TFile(offsetFileName.data(),"read");
  }  else { 
    Recenter_ = false;
  }
  mx = ntrkbins;
  if(!useNtrk_) {
    mx = ncentbins;
  }
  
  std::cout<<"==============================================="<<std::endl;
  if(genMC_) std::cout<<"Generator MC Replay"<<std::endl;
  std::cout<<"centralityBinTag_           "<<centralityBinTag_.encode()<<std::endl;
  std::cout<<"centralityTag_              "<<centralityTag_.encode()<<std::endl;
  std::cout<<"vertexTag_                  "<<vertexTag_.encode()<<std::endl;
  std::cout<<"trackTag_                   "<<trackTag_.encode()<<std::endl;
  std::cout<<"inputPlanesTag_             "<<inputPlanesTag_.encode()<<std::endl;
  std::cout<<"FlatOrder_                  "<<FlatOrder_<<std::endl;
  std::cout<<"NumFlatBins_                "<<NumFlatBins_<<std::endl;
  std::cout<<"flatnvtxbins_               "<<flatnvtxbins_<<std::endl;
  std::cout<<"flatminvtx_                 "<<flatminvtx_<<std::endl;
  std::cout<<"flatdelvtx_                 "<<flatdelvtx_<<std::endl;
  std::cout<<"caloCentRef_                "<<caloCentRef_<<std::endl;
  std::cout<<"caloCentRefWidth_           "<<caloCentRefWidth_<<std::endl;
  std::cout<<"CentBinCompression_         "<<CentBinCompression_<<std::endl;
  std::cout<<"Noffmin_                    "<<Noffmin_<<std::endl;
  std::cout<<"Noffmax_                    "<<Noffmax_<<std::endl;
  std::cout<<"minrun_                     "<<minrun_<<std::endl;
  std::cout<<"maxrun_                     "<<maxrun_<<std::endl;
  std::cout<<"effTable_                   "<<effTable_<<std::endl;
  std::cout<<"dzerror_                    "<<dzdzerror_<<endl;
  std::cout<<"d0error_                    "<<d0d0error_<<endl;
  std::cout<<"pterror_                    "<<pterror_<<endl;
  std::cout<<"nvtx_                       "<<nvtx_<<endl;
  if(bCaloMatching_) { 
    std::cout<<"bCaloMatching_              true"<<std::endl;
    std::cout<<"reso_                     "<<reso_<<std::endl;   
  }
  std::cout<<"dzdzerror_pix_               "<<dzdzerror_pix_<<std::endl;
  std::cout<<"chi2_                        "<<chi2_<<std::endl;
  if(Recenter_) { 
    std::cout<<"Recenter_                  true"<<std::endl;
  } else {
    std::cout<<"Recenter_                  false"<<std::endl;
  }
  if(useNtrk_) { 
    std::cout<<"Ntrk binning                  "<<std::endl;
  } else {
    std::cout<<"Centrality binning            "<<std::endl;
  }
  std::cout<<"==============================================="<<std::endl;
  TDirectory * save = gDirectory;
  TFileDirectory conddir = fs->mkdir("Conditions");
  conddir.make<TH1I>(centralityBinTag_.label().data(),centralityBinTag_.label().data(),1,0,1);
  conddir.make<TH1I>(centralityTag_.label().data(), centralityTag_.label().data(),1,0,1);
  conddir.make<TH1I>(vertexTag_.label().data(), vertexTag_.label().data(),1,0,1);
  conddir.make<TH1I>(trackTag_.label().data(), trackTag_.label().data(),1,0,1);
  conddir.make<TH1I>(inputPlanesTag_.label().data(), inputPlanesTag_.label().data(),1,0,1);
  string etable = Form("EffTable_%s",effTable_.data());
  conddir.make<TH1I>(etable.data(), etable.data(),1,0,1);
  string note_FlatOrder = Form("FlatOrder_%d",FlatOrder_);
  conddir.make<TH1I>(note_FlatOrder.data(), note_FlatOrder.data(),1,0,1);
  string note_NumFlatBins = Form("NumFlatBins_%d",NumFlatBins_);
  conddir.make<TH1I>(note_NumFlatBins.data(), note_NumFlatBins.data(),1,0,1);
  string note_flatnvtxbins = Form("flatnvtxbins_%d",flatnvtxbins_);
  conddir.make<TH1I>(note_flatnvtxbins.data(), note_flatnvtxbins.data(),1,0,1);
  string note_flatminvtx = Form("flatminvtx_%07.2f",flatminvtx_);
  conddir.make<TH1I>(note_flatminvtx.data(), note_flatminvtx.data(),1,0,1);
  string note_flatdelvtx = Form("flatdelvtx_%07.2f",flatdelvtx_);
  conddir.make<TH1I>(note_flatdelvtx.data(), note_flatdelvtx.data(),1,0,1);
  string note_caloCentRef = Form("caloCentRef_%d",(int)caloCentRef_);
  conddir.make<TH1I>(note_caloCentRef.data(), note_caloCentRef.data(),1,0,1);
  string note_caloCentRefWidth = Form("caloCentRefWidth_%d",(int)caloCentRefWidth_);
  conddir.make<TH1I>(note_caloCentRefWidth.data(), note_caloCentRefWidth.data(),1,0,1);
  string note_dzdzerror = Form("dzdzerror_%07.2f",dzdzerror_);
  conddir.make<TH1I>(note_dzdzerror.data(), note_dzdzerror.data(),1,0,1);
  string note_d0d0error = Form("d0d0error_%07.2f",d0d0error_);
  conddir.make<TH1I>(note_d0d0error.data(), note_d0d0error.data(),1,0,1);
  string note_dzdzerror_pix = Form("dzdzerror_pix_%07.2f",dzdzerror_pix_);
  conddir.make<TH1I>(note_dzdzerror_pix.data(), note_dzdzerror_pix.data(),1,0,1);
  string note_chi2 = Form("chi2_%07.2f",chi2_);
  conddir.make<TH1I>(note_chi2.data(), note_chi2.data(),1,0,1);
  string note_vtx_range = Form("vtx_%5.1f_%5.1f",minvz_,maxvz_);
  conddir.make<TH1I>(note_vtx_range.data(), note_vtx_range.data(),1,0,1);
  string note_nvtx = Form("nvtx_%d",nvtx_);
  conddir.make<TH1I>(note_nvtx.data(), note_nvtx.data(),1,0,1);
  if(bCaloMatching_) {
    conddir.make<TH1I>("bCaloMatching_Set_True", "bCaloMatching_Set_True",1,0,1);
    string note_reso = Form("reso_%07.2f",reso_);
    conddir.make<TH1I>(note_reso.data(),note_reso.data(),1,0,1);
  } else {
    conddir.make<TH1I>("bCaloMatching_Set_False", "bCaloMatching_Set_False",1,0,1);
  }
  if(Recenter_) {
    conddir.make<TH1I>("RecenterTracks", "RecenterTracks",1,0,1);
  } else {
    conddir.make<TH1I>("DoNotRecenterTracks", "DoNotRecenterTracks",1,0,1);
  }
  if(useNtrk_) {
    conddir.make<TH1I>("useNtrk_Set_True", "useNtrk_Set_True",1,0,1);
  } else {
    conddir.make<TH1I>("useNtrk_Set_False", "useNtrk_Set_False",1,0,1);
  }
  
  save->cd();
  hNtrk = fs->make<TH1D>("Ntrk","Ntrk",1001,0,3000);
  hNoff = fs->make<TH1D>("Noff","Noff",1001,0,3000);
  int npt = nptbins;
  for(int iorder = 1; iorder<=7; iorder++) {
    qxtrk[iorder-1] = fs->make<TH2F>(Form("qxtrk%d",iorder),Form("qxtrk%d",iorder),npt,ptbins, netabinsDefault, etabinsDefault);
    qytrk[iorder-1] = fs->make<TH2F>(Form("qytrk%d",iorder),Form("qytrk%d",iorder),npt,ptbins, netabinsDefault, etabinsDefault);
    qxtrk[iorder-1]->SetOption("colz");
    qytrk[iorder-1]->SetOption("colz");
    qxtrk[iorder-1]->Sumw2();
    qytrk[iorder-1]->Sumw2();
    qxtrk[iorder-1]->SetXTitle("p_{T} (GeV/c");
    qxtrk[iorder-1]->SetYTitle(Form("#eta (n=%d)",iorder));
    qytrk[iorder-1]->SetXTitle("p_{T} (GeV/c");
    qytrk[iorder-1]->SetYTitle(Form("#eta (n=%d)",iorder));
  }
  qcnt =  fs->make<TH2F>("qcnt", "qcnt",npt,ptbins, netabinsDefault, etabinsDefault);
  qcnt->SetXTitle("p_{T} (GeV/c");
  qcnt->SetYTitle("#eta");
  avpt =  fs->make<TH2F>("avpt","avpt",npt,ptbins, netabinsDefault, etabinsDefault);
  qcnt->SetOption("colz");
  avpt->SetOption("colz");
  qcnt->Sumw2();
  avpt->Sumw2();
  hTemplate = (TH2D *) qcnt->Clone("hTemplate");
  hTemplate->SetDirectory(0);
  hTemplate->Reset();
  if(!trackq->is_GenMC()) {
    for(int i = 0; i<mx; i++) {
      for(int j = 1; j<=7; j++){
	for(int k = 0; k<flatnvtxbins_; k++) {
	  offsets[k].wqxtrkRef[j-1][i] = (TH2D *) frecenter->Get(Form("wqxtrk_%d_%d_%d",j,k,i));
	  offsets[k].wqytrkRef[j-1][i] = (TH2D *) frecenter->Get(Form("wqytrk_%d_%d_%d",j,k,i));
	  offsets[k].wqxtrkRef[j-1][i]->SetDirectory(0);
	  offsets[k].wqytrkRef[j-1][i]->SetDirectory(0);
	}
      }
    }
  } else {
    for(int i = 0; i<mx; i++) {
      for(int j = 1; j<=7; j++){
	for(int k = 0; k<flatnvtxbins_; k++) {
	  offsets[k].wqxtrkRef[j-1][i] = (TH2D *) hTemplate->Clone(Form("wqxtrk_%d_%d_%d",j,k,i));
	  offsets[k].wqytrkRef[j-1][i] = (TH2D *) hTemplate->Clone(Form("wqytrk_%d_%d_%d",j,k,i));
	  offsets[k].wqxtrkRef[j-1][i]->Reset();
	  offsets[k].wqytrkRef[j-1][i]->Reset();
	  offsets[k].wqxtrkRef[j-1][i]->SetDirectory(0);
	  offsets[k].wqytrkRef[j-1][i]->SetDirectory(0);
	}
      }
    }
  }
  hcent = fs->make<TH1D>("cent","cent",220,-10,110);
  hvtx = fs->make<TH1D>("vtx","vtx",600,-30,30);
  hvtxRaw = fs->make<TH1D>("vtxRaw","vtxRaw",600,-30,30);
  hb = fs->make<TH1D>("b","b",10000,0,10000);
  hRP = fs->make<TH1D>("rp","rp",1000,-4,4);
  hcentbins = fs->make<TH1D>("centbins","centbins",201,0,200);
  if(useNtrk_) {
    hcentres = fs->make<TH1D>("centres","centres",ntrkbins,trkBins);
  } else {
    hcentres = fs->make<TH1D>("centres","centres",ncentbins,centbins);
  }
  hrun = fs->make<TH1I>("runs","runs",maxrun_-minrun_+1,minrun_,maxrun_);
  hptNtrk = fs->make<TH1D>("ptNtrk","ptNtrk",npt,ptbins);
  hptNtrk->SetXTitle("p_{T} (GeV/c)");
  hptNtrk->SetYTitle("Ntrks (|#eta|<1; 0-5)");
  hptNtrkGood = fs->make<TH1D>("ptNtrkGood","ptNtrkGood",npt,ptbins);
  hptNtrkGood->SetXTitle("p_{T} (GeV/c)");
  hptNtrkGood->SetYTitle("Ntrks (Good) (|#eta|<1; 0-5)");
  hNtrkRet = fs->make<TH1I>("NtrkRet","NtrkRet", 10,0,10);
  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";
  NumFlatBins_ = ntrkbins;
  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    TFileDirectory subdir = fs->mkdir(Form("%s",EPNames[i].data()));
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->init(FlatOrder_,NumFlatBins_,flatnvtxbins_,flatminvtx_,flatdelvtx_,EPNames[i],EPOrder[i]);
    Double_t psirange = 4;
    if(EPOrder[i]==1 ) psirange = 3.5;
    if(EPOrder[i]==2 ) psirange = 2;
    if(EPOrder[i]==3 ) psirange = 1.5;
    if(EPOrder[i]==4 ) psirange = 1;
    if(EPOrder[i]==5) psirange = 0.8;
    if(EPOrder[i]==6) psirange = 0.6;
    if(EPOrder[i]==7) psirange = 0.5;
    
    hPsi[i] = subdir.make<TH1D>("psi","psi",800,-psirange,psirange);
    hPsi[i]->SetXTitle("#Psi");
    hPsi[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
    hPsiOffset[i] = subdir.make<TH1D>("psiOffset","psiOffset",800,-psirange,psirange);
    hPsiOffset[i]->SetXTitle("#Psi");
    hPsiOffset[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
    
    hPsiFlat[i] = subdir.make<TH1D>("psiFlat","psiFlat",800,-psirange,psirange);
    hPsiFlat[i]->SetXTitle("#Psi");
    hPsiFlat[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
  }
  int nanalbins = 0;
  if(useNtrk_) {
    nanalbins = ntrkbins;
  } else {
    nanalbins = ncentbins;
  }
  //============= MC =========================
  if(genMC_) {
    TFileDirectory mcdir = fs->mkdir("MC");
    for(int i = 0; i<nanalbins; i++) {
      TFileDirectory mcsubdir; 
      if(useNtrk_) {
	mcsubdir = mcdir.mkdir(Form("%d_%d",(int)trkBins[i],(int)trkBins[i+1]));
      } else {
	mcsubdir = mcdir.mkdir(Form("%d_%d",(int)centbins[i],(int)centbins[i+1]));
      }
    }
  }
  //==============   Resolution terms  ========
  TFileDirectory resdir = fs->mkdir("Resolutions");
  for(int i = 0; i<nanalbins; i++) {
    TFileDirectory ressubdir;
    if(useNtrk_) {
      ressubdir = resdir.mkdir(Form("%d_%d",(int)trkBins[i],(int)trkBins[i+1]));
    } else {
      ressubdir = resdir.mkdir(Form("%d_%d",(int)centbins[i],(int)centbins[i+1]));
    }
    for(int iorder = 1; iorder<=7; iorder++) {
      res[iorder-1][i] = ressubdir.make<TH2D>(Form("res%d",iorder),Form("res%d",iorder),50,0,50,50,0,50);
      resw[iorder-1][i] = ressubdir.make<TH2D>(Form("resw%d",iorder),Form("resw%d",iorder),50,0,50,50,0,50);
      resep[iorder-1][i] = ressubdir.make<TH2D>(Form("resep%d",iorder),Form("resep%d",iorder),50,0,50,50,0,50);
      rescnt[iorder-1][i] = ressubdir.make<TH2D>(Form("rescnt%d",iorder),Form("rescnt%d",iorder),50,0,50,50,0,50);
      res[iorder-1][i]->Reset();
      res[iorder-1][i]->Sumw2();
      res[iorder-1][i]->SetOption("colz");
      resw[iorder-1][i]->Reset();
      resw[iorder-1][i]->Sumw2();
      resw[iorder-1][i]->SetOption("colz");
      resep[iorder-1][i]->Reset();
      resep[iorder-1][i]->Sumw2();
      resep[iorder-1][i]->SetOption("colz");
      rescnt[iorder-1][i]->Reset();
      rescnt[iorder-1][i]->Sumw2();
      rescnt[iorder-1][i]->SetOption("colz");
    }
    res12[i] = ressubdir.make<TH2D>("res12","res12",50,0,50,50,0,50);
    resw12[i] = ressubdir.make<TH2D>("resw12","resw12",50,0,50,50,0,50);
    resep12[i] = ressubdir.make<TH2D>("resep12","resep12",50,0,50,50,0,50);
    rescnt12[i] = ressubdir.make<TH2D>("rescnt12","rescnt12",50,0,50,50,0,50);
    res12[i]->Reset();
    resw12[i]->Reset();
    resep12[i]->Reset();
    rescnt12[i]->Reset();
    res12[i]->Sumw2();
    resw12[i]->Sumw2();
    resep12[i]->Sumw2();
    rescnt12[i]->Sumw2();
    res12[i]->SetOption("colz");
    resw12[i]->SetOption("colz");
    resep12[i]->SetOption("colz");
    rescnt12[i]->SetOption("colz");
    res12[i]->GetXaxis()->SetTitle("n = 1");
    res12[i]->GetYaxis()->SetTitle("n = 2");
    resw12[i]->GetXaxis()->SetTitle("n = 1");
    resw12[i]->GetYaxis()->SetTitle("n = 2");
    resep12[i]->GetXaxis()->SetTitle("n = 1");
    resep12[i]->GetYaxis()->SetTitle("n = 2");
    rescnt12[i]->GetXaxis()->SetTitle("n = 1");
    rescnt12[i]->GetYaxis()->SetTitle("n = 2");
  }
  //====================================================
  //========= Spectra =================================
  TFileDirectory specdir = fs->mkdir("Spectra");
  for(int i = 0; i<nanalbins; i++) {
    TFileDirectory specsubdir;
    if(useNtrk_) {
      specsubdir = specdir.mkdir(Form("%d_%d",(int)trkBins[i],(int)trkBins[i+1]));
    } else {
      specsubdir = specdir.mkdir(Form("%d_%d",(int)centbins[i],(int)centbins[i+1]));
    }
    ptspec[i] =  specsubdir.make<TH2D>("ptspec","ptspec",npt,ptbins, netabinsDefault, etabinsDefault);
    ptspec[i]->Sumw2();
    ptspec[i]->SetXTitle("p_{T} (GeV/c)");
    ptspec[i]->SetYTitle("Yield");
    ptspec[i]->SetOption("colz");
  }
  
  //==============================
  if(makeTree_) {
    tree = fs->make<TTree>("tree","EP tree");
    tree->Branch("Cent",&centval,"cent/D");
    tree->Branch("NtrkOff",&Noff,"Noff/I");
    tree->Branch("ntrkflat",&ntrkval,"nofftrak/I");
    tree->Branch("Vtx",&vtx,"vtx/D");
    tree->Branch("epang",&epang, epnames.Data());
    tree->Branch("eporig",&eporig, epnames.Data());
    tree->Branch("qx",      &qx,       epnames.Data());
    tree->Branch("qy",      &qy,       epnames.Data());
    tree->Branch("q",       &q,       epnames.Data());
    tree->Branch("vn", &vn, epnames.Data());
    tree->Branch("mult",    &epmult,  epnames.Data());
    tree->Branch("sumw",    &sumw,  epnames.Data());
    tree->Branch("sumw2",    &sumw2,  epnames.Data());
    tree->Branch("Run",     &runno_,   "run/i");
    tree->Branch("Rescor",  &rescor,   epnames.Data());
    tree->Branch("RescorErr",  &rescorErr,   epnames.Data());
    tree->Branch("qxtrk1",  "TH2F",  &qxtrk[0], 128000, 0);
    tree->Branch("qytrk1",  "TH2F",  &qytrk[0], 128000, 0);
    tree->Branch("qxtrk2",  "TH2F",  &qxtrk[1], 128000, 0);
    tree->Branch("qytrk2",  "TH2F",  &qytrk[1], 128000, 0);
    tree->Branch("qxtrk3",  "TH2F",  &qxtrk[2], 128000, 0);
    tree->Branch("qytrk3",  "TH2F",  &qytrk[2], 128000, 0);
    tree->Branch("qxtrk4",  "TH2F",  &qxtrk[3], 128000, 0);
    tree->Branch("qytrk4",  "TH2F",  &qytrk[3], 128000, 0);
    tree->Branch("qxtrk5",  "TH2F",  &qxtrk[4], 128000, 0);
    tree->Branch("qytrk5",  "TH2F",  &qytrk[4], 128000, 0);
    tree->Branch("qxtrk6",  "TH2F",  &qxtrk[5], 128000, 0);
    tree->Branch("qytrk6",  "TH2F",  &qytrk[5], 128000, 0);
    tree->Branch("qxtrk7",  "TH2F",  &qxtrk[6], 128000, 0);
    tree->Branch("qytrk7",  "TH2F",  &qytrk[6], 128000, 0);
    tree->Branch("qcnt",    "TH2F",  &qcnt, 128000, 0);
    tree->Branch("avpt",    "TH2F",  &avpt, 128000, 0);
  }
}


VNAnalyzer::~VNAnalyzer()
{
  if(frecenter!=NULL) frecenter->Close();  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VNAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  Bool_t newrun = kFALSE;
  if(runno_ != iEvent.id().run()) newrun = kTRUE;
  runno_ = iEvent.id().run();
  hrun->Fill(runno_);
  if(FirstEvent_ || newrun) {
    FirstEvent_ = kFALSE;
    newrun = kFALSE;
    //
    //Get Size of Centrality Table
    //
    nCentBins_ = 200;
    if(!useNtrk_ && !genMC_) {
      edm::ESHandle<CentralityTable> centDB_;
      iSetup.get<HeavyIonRcd>().get(centralityLabel_,centDB_);
      nCentBins_ = (int) centDB_->m_table.size();
      for(int i = 0; i<NumEPNames; i++) {
    	flat[i]->setCaloCentRefBins(-1,-1);
    	if(caloCentRef_>0) {
    	  int minbin = (caloCentRef_-caloCentRefWidth_/2.)*nCentBins_/100.;
    	  int maxbin = (caloCentRef_+caloCentRefWidth_/2.)*nCentBins_/100.;
    	  minbin/=CentBinCompression_;
    	  maxbin/=CentBinCompression_;
    	  if(minbin>0 && maxbin>=minbin) {
    	    if(EPDet[i]==HF || EPDet[i]==Castor) flat[i]->setCaloCentRefBins(minbin,maxbin);
    	  }
    	}
      }
    }
    //
    // Get flattening parameter file.  
    //
    if(!genMC_) {
      edm::ESHandle<RPFlatParams> flatparmsDB_;
      iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
      LoadEPDB * db = new LoadEPDB(flatparmsDB_,flat);
      if(!db->IsSuccess()) {
	std::cout<<"Flattening db inconsistancy, will continue without: "<<std::endl;
	loadDB_ = kFALSE;
      }        
    }
  } //First event
  
  //
  // Get Centrality
  //

  int Noff=0;
 
  int bin = 0;
  int cbin = 0;
  if(!useNtrk_) {
    Noff = getNoff(iEvent,iSetup);
    hNoff->Fill(Noff);
    iEvent.getByToken(centralityBinToken, cbin_);
    cbin = *cbin_;
    bin = cbin/CentBinCompression_; 
    double cscale = 100./nCentBins_;
    centval = cscale*cbin;
    hcentres->Fill(centval);
  } else {
    
    iEvent.getByToken(centralityToken, centrality_);
    //ntrkval = centrality_->Ntracks();
    Noff = getNoff(iEvent, iSetup);
    hNoff->Fill(Noff);
    bin = NtrkToBin(Noff)-1;
    cbin = bin;
    centval = Noff;
    hcentres->Fill(Noff);
  }
  int ibin = hcentres->FindBin(centval)-1;
  hcent->Fill(centval);
  hcentbins->Fill(cbin);
  if(genMC_) {
    iEvent.getByToken(bToken, b_);
    b = *b_;     
    hb->Fill(500.*b);
    iEvent.getByToken(RPToken, RP_);
    rp = *RP_;     
    hRP->Fill(rp);

  } 
  int ntrkval=fillTracks(iEvent, iSetup, ibin);
  if(ntrkval<=0) return;
  hvtx->Fill(vtx);
  hNtrk->Fill(ntrkval);
  //
  // Get Event Planes
  //
  iEvent.getByToken(inputPlanesToken,inputPlanes_);
  
  if(!inputPlanes_.isValid()){
    cout << "Error! Can't get hiEvtPlaneFlat product!" << endl;
    return ;
  }
  
  Int_t indx = 0;
  for(int i = 0; i<NumEPNames; i++) {
    epang[i] = -10;
    epsin[i] = 0;
    epcos[i] = 0;
    qx[i] = 0;
    qy[i] = 0;
    q[i] = 0;
    vn[i] = 0;
    epmult[i] = 0;
    sumw[i] = 0;
    sumw2[i] = 0;
  }
  for (EvtPlaneCollection::const_iterator rp = inputPlanes_->begin();rp !=inputPlanes_->end(); rp++) {
    if(indx != rp->indx() ) {
      cout<<"EP inconsistency found. Abort."<<endl;
      return;
    }
    if(rp->sumSin()!=0 || rp->sumCos()!=0) {
      if(rp->mult()>3 && fabs(vtx)<15) {
	epang[indx]=rp->angle();
	eporig[indx]=rp->angle(0);
	epsin[indx] = rp->sumSin();
	epcos[indx] = rp->sumCos();
	hPsi[indx]->Fill(rp->angle(0));
	hPsiOffset[indx]->Fill(rp->angle(1));
	hPsiFlat[indx]->Fill(rp->angle(2));
	qx[indx]  = rp->qx(); 
	qy[indx]  = rp->qy();
	q[indx]   = rp->q();
	vn[indx] = rp->vn(0);
	epmult[indx] = (double) rp->mult();
	sumw[indx] = rp->sumw();
	sumw2[indx] = rp->sumw2();
	if(!useNtrk_) {
	  rescor[indx] = flat[indx]->getCentRes1((int) centval);
	  rescorErr[indx] = flat[indx]->getCentResErr1((int) centval);
	} else {
	  rescor[indx] = 0.;
	  rescorErr[indx] = 0.;
	}
      } 
    }
    ++indx; 
  }
  
  for(int iorder = 1; iorder <=7; iorder++) {
    int epmin = 0;
    int epmax = 0;
    if(iorder == 1 ) {
      epmin = HFm1;
      epmax = HFp1f;
    }else if (iorder == 2) {
      epmin = HFm2;
      epmax = HFp2f;
    }else if (iorder == 3) {
      epmin = HFm3;
      epmax = HFp3f;
    }else if (iorder == 4) {
      epmin = HFm4;
      epmax = HFp4f;
    }else if (iorder == 5) {
      epmin = HFm5;
      epmax = trackp522;
    }else if (iorder == 6) {
      epmin = HFm6;
      epmax = trackp622;
    }else if (iorder == 7) {
      epmin = HFm7;
      epmax = trackp722;
    }
    for(int i = epmin; i<= epmax; i++) {
      for(int j = i; j<=epmax; j++) {
	double ii = i-epmin;
	double jj = j-epmin;
	double w =(pow(qx[i],2)+pow(qy[i],2))*(pow(qx[j],2)+pow(qy[j],2));
	if(w>0) {
	  res[iorder-1][ibin]->Fill(ii,jj,   qx[i]*qx[j]+qy[i]*qy[j]);
	  resw[iorder-1][ibin]->Fill(ii,jj,  sumw[i]*sumw[j]);
	  resep[iorder-1][ibin]->Fill(ii,jj, (qx[i]*qx[j]+qy[i]*qy[j])/sqrt(w));
	  rescnt[iorder-1][ibin]->Fill(ii,jj);
	}
      }
    }
    for(int i = HFm1; i<= HFp1f; i++) {
      for(int j = HFm2; j<=HFp2f; j++) {
	double ii = i-HFm1;
	double jj = j-HFm2;
	double w =pow((pow(qx[i],2)+pow(qy[i],2)),2)*(pow(qx[j],2)+pow(qy[j],2));
	if(w>0) {
	  res12[ibin]->Fill(ii,jj,  pow(qx[i],2)*qx[j]-pow(qy[i],2)*qy[j]+2*qx[i]*qy[i]*qy[j]);
	  resw12[ibin]->Fill(ii,jj,  pow(sumw[i],2)*sumw[j]);
	  resep12[ibin]->Fill(ii,jj,  ( pow(qx[i],2)*qx[j]-pow(qy[i],2)*qy[j]+2*qx[i]*qy[i]*qy[j])/sqrt(w));
	  rescnt12[ibin]->Fill(ii,jj);
	}
      }
    }

  }
  if(makeTree_) tree->Fill(); 
}



// ------------ method called once each job just before starting event loop  ------------
void 
VNAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VNAnalyzer::endJob() {
}
bool
VNAnalyzer::CaloMatch(const reco::Track & track, const edm::Event & iEvent, unsigned int idx)
{
  if ( !bCaloMatching_ ) return true;
  edm::Handle<reco::PFCandidateCollection> pfCand;
  iEvent.getByToken( pfToken_, pfCand );
  double energy = 0;
  for ( reco::PFCandidateCollection::const_iterator it = pfCand->begin(); it != pfCand->end(); ++it ) {
    reco::TrackRef trackRef = it->trackRef();
    if ( !((it->particleId() != reco::PFCandidate::h) ||
	   (it->particleId() != reco::PFCandidate::e) ||
	   (it->particleId() != reco::PFCandidate::mu) )) continue;
    if ( idx == trackRef.key() ) {
      energy = it->ecalEnergy() + it->hcalEnergy();
      break;
    }
  }
  
  if( track.pt() < 20 || ( energy/( track.pt()*TMath::CosH(track.eta() ) ) > reso_ && (energy)/(TMath::CosH(track.eta())) > (track.pt() - 80.0) )  ) return true;
  else {
    return false;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(VNAnalyzer);

