import FWCore.ParameterSet.Config as cms

evtPlaneCalibTree = cms.EDAnalyzer("EvtPlaneCalibTree",
                            vertexTag_=cms.InputTag("offlinePrimaryVerticesRecovery"),
                            centralityTag_=cms.InputTag("hiCentrality"),
                            caloTag_ = cms.InputTag("towerMaker"),
                            centralityBinTag_ = cms.InputTag("centralityBin","HFtowers"),
                            centralityVariable = cms.string("HFtowers"),
                            inputPlanesTag_ = cms.InputTag("hiEvtPlane",""),
                            FlatOrder_ = cms.untracked.int32(9),
                            NumFlatBins_ = cms.untracked.int32(40),
                            CentBinCompression_ = cms.untracked.int32(5),
                            minet_ = cms.untracked.double(0.01),
                            maxet_ = cms.untracked.double(-1.),
                            minpt_ = cms.untracked.double(0.5),
                            maxpt_ = cms.untracked.double(3.0),
                            flatnvtxbins_ = cms.int32(10),
                            flatminvtx_ = cms.double(-15.0),
                            flatdelvtx_ = cms.double(3.0),
                            chi2perlayer_ = cms.untracked.double(0.18),
                            dzdzerror_ = cms.untracked.double(3.),
                            d0d0error_ = cms.untracked.double(3.),
                            pterror_ = cms.untracked.double(0.1),
                            useNtrkBins_ = cms.untracked.bool(False),
                            genMC_ = cms.untracked.bool(False),
                            bTag_ = cms.InputTag("mcEvtPlane","b","FlatCalib"),
                            bypassCentrality_ = cms.untracked.bool(False),
                            trackTag = cms.InputTag("generalTracks"),
                            minvz_ = cms.untracked.double(-15.),
                            maxvz_ = cms.untracked.double(15.),
                            dzdzerror_pix_ = cms.untracked.double(40.0),
                            chi2_ = cms.untracked.double(40.)
                            )
                            




    
