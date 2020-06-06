import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
ivars = VarParsing.VarParsing('standard')

ivars.register ('lumifile',
                'Cert_326381-327560_HI_PromptReco_Collisions18_JSON.txt',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="lumi file")

ivars.register ('offset',
                'offset_PbPb2018_merged_1_326544.root',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="offset file")

ivars.register ('dbfile',
                'HeavyIonRPRcd_PbPb2018_March2020_offline.db',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="dbfile file")

ivars.register ('eff',                'NULL',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="efficiency file")

ivars.register('tracks',
                'generalTracks',
		VarParsing.VarParsing.multiplicity.singleton,
		VarParsing.VarParsing.varType.string,
                "track collection")

ivars.register('part',
                'KS',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "CH, LM or KS")

ivars.register('massRange',
                'Peak',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "Peak/SBPos/SBNeg")

ivars.register('rap',
                'Mid',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "Mid/Fwd")

ivars.register('BDT',
                0.15,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.float,
                "BDT cut")


ivars.parseArguments()

massLow = 0.;
massHigh = 0.;

if ivars.part == 'LM':
    if ivars.massRange == 'Peak':
        massLow = 1.1115
        massHigh = 1.1200
    if ivars.massRange == 'SBPos':
        massLow = 1.122
        massHigh = 1.16
    if ivars.massRange == 'SBNeg':
        massLow = 1.08
        massHigh = 1.1095

if ivars.part == 'KS':
    if ivars.massRange == 'Peak':
        massLow = 0.492
        massHigh = 0.503
    if ivars.massRange == 'SBPos':
        massLow = 0.506
        massHigh = 0.565
    if ivars.massRange == 'SBNeg':
        massLow = 0.43
        massHigh = 0.489

print ivars.part
print ivars.massRange
print ivars.rap
print ivars.BDT

process = cms.Process("VNANAL")
process.load('Configuration.StandardSequences.Services_cff')
process.load("CondCore.CondDB.CondDB_cfi")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")
process.load('HeavyIonsAnalysis.EventAnalysis.clusterCompatibilityFilter_cfi')
process.load("HeavyIonsAnalysis.Configuration.hfCoincFilter_cff")
process.load("HeavyIonsAnalysis.Configuration.analysisFilters_cff")
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")


process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")
process.load("HeavyIonsAnalysis.VNAnalysis/vnanalyzer_cfi")


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_v2', '')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1033p1x01_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")

process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('RecoHI.HiCentralityAlgos.CentralityFilter_cfi')
process.load('PbPb18_HIMB_rereco')

process.eventSelection = cms.Sequence(
	process.primaryVertexFilter
	* process.hfCoincFilter2Th4
	* process.clusterCompatibilityFilter
    )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30000) )
process.MessageLogger.cerr.FwkReport.reportEvery=10000

process.CondDB.connect = "sqlite_file:"+ivars.dbfile
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                       process.CondDB,
                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
#                                                                  tag = cms.string('HeavyIonRPRcd_PbPb2018_offline')
                                                                  tag = cms.string('HeavyIonRPRcd')
                                                                  )
                                                         )
                                      )
process.es_prefer_flatparms = cms.ESPrefer('PoolDBESSource','')


import FWCore.PythonUtilities.LumiList as LumiList
goodLumiSecs = LumiList.LumiList(filename = ivars.lumifile ).getCMSSWString().split(',')

if ivars.part == 'CH':
    print 'CH data stream'
    process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F873AB58-6232-974D-9CDB-3BFB1504B449.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F7C4BAFB-F443-6C4C-8FFF-B2C89EA3A776.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F6200B99-E688-8E47-8301-6FB2D0521021.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F40A87C5-494C-0248-BC2F-AB9442E1DDD6.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F3CF4F88-F926-224E-9B77-3BF0B3F084A6.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F31298D8-F6D5-E941-B9A0-DD76272B9AE0.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F23CAED0-68CB-EC46-8EEB-1EAACC9556A6.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F23AC73A-E593-6B4E-BA9D-15D55D6EC5FC.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F1C264FD-6BF9-294B-8619-8277686EA99A.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F0B44001-C3AA-D640-AAF1-B3417A46D83D.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/F04765AE-CF2A-9643-909D-A108C0025854.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/EFBF111B-B428-FB44-B620-615A55E263C4.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/EFA366C6-3F15-0149-B0F2-D456885C9BA7.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/EE2C2C5B-43D5-BA44-A465-743A22570095.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/ED7B7372-79AB-BC4F-83A0-E5027732081E.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/ED1D69A8-BD38-F648-AA4F-B15F95D5EBCD.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/ED0E1364-EFA6-154D-B9C5-4F8A171926E3.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/ECAA8568-3828-E24B-97DD-886F9B800DFD.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/EC5143AB-0831-FA4C-8919-3DDA127D4258.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/EA3951F7-FE50-7545-947A-9577298DC61B.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/EA129831-A476-2841-888C-3B1ECDDEB7FB.root',
            'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias1/AOD/04Apr2019-v1/610007/E7B62111-338A-714F-B540-CED49C0C017F.root'),
                                 inputCommands=cms.untracked.vstring(
            'keep *',
            'drop *_hiEvtPlane_*_*'
            )
                                 )
else :
    print 'v0 data stream'
    process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring("file:PbPb2018_RERECO_V0Skim.root"),
                                secondaryFileNames = cms.untracked.vstring(
            "file:FF31F840-542E-1A49-ACF7-9043F8169E67.root"
            ),
                                )


import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltMB = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMB.HLTPaths = [
	"HLT_HIMinimumBias_*"
	]
process.hltMB.andOr = cms.bool(True)
process.hltMB.throw = cms.bool(False)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("vnanal.root")
)


if ivars.part == 'LM' or ivars.part == 'KS' :

    if ivars.part == 'LM':
        process.QWV0Event = process.QWV0EventLm.clone()
    if ivars.part == 'KS':
            process.QWV0Event = process.QWV0EventKs.clone()
            
    if ivars.rap == 'Fwd':
        process.QWV0Event.cuts[0].AbsRapmax = cms.untracked.double(2.0)
        process.QWV0Event.cuts[0].AbsRapmin = cms.untracked.double(1.0)
        
        
    process.QWV0MVA = cms.EDProducer('QWV0MVAVectProducer'
                                     , vertexSrc = cms.untracked.InputTag('offlinePrimaryVerticesRecovery')
                                     , V0Src = cms.untracked.InputTag('QWV0Event', 'Lambda')
                                     , dbCent = cms.untracked.InputTag('dbCent')
                                     , mvaXML = cms.untracked.string('MC_Full_BDT_D4.LM.weights.xml')
                                     , mvaCut = cms.untracked.double(ivars.BDT)
                                     )
        
        
    if ivars.part == 'KS':
        process.QWV0MVA.V0Src = cms.untracked.InputTag('QWV0Event', 'Kshort')
        process.QWV0MVA.mvaXML = cms.untracked.string('MC_Full_BDT_D4.KS.weights.xml')
        
    process.QWV0MVAVector = cms.EDProducer('QWV0VectProducer'
                                           , vertexSrc = cms.untracked.InputTag('offlinePrimaryVerticesRecovery')
                                           , trackSrc = cms.untracked.InputTag('generalTracks')
                                           , V0Src = cms.untracked.InputTag('QWV0MVA', 'Lambda')
                                           , daughter_cuts = cms.untracked.PSet(
            )
                                           , cuts = cms.untracked.VPSet(
            cms.untracked.PSet(
                )
            )
                                           )
    
    if ivars.part == 'KS':
        process.QWV0MVAVector.V0Src = cms.untracked.InputTag('QWV0MVA', 'Kshort')
        
        
process.QWEvent.Year = cms.untracked.int32(2018)
process.QWEvent.chi2  = cms.untracked.double(0.18)
        
process.dump = cms.EDAnalyzer("EventContentAnalyzer")
        
process.hiEvtPlane.trackTag = cms.InputTag(ivars.tracks)
process.hiEvtPlane.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlane.loadDB = cms.bool(True)
process.hiEvtPlane.useNtrk = cms.untracked.bool(False)


process.hiEvtPlaneFlat.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlaneFlat.useNtrk = cms.untracked.bool(False)
process.vnanalyzer.trackTag_ = cms.InputTag(ivars.tracks)
process.vnanalyzer.vertexTag_ = cms.InputTag("offlinePrimaryVerticesRecovery")
process.vnanalyzer.useNtrk = cms.untracked.bool(False)
process.vnanalyzer.offsetFile = cms.untracked.string( ivars.offset )
process.vnanalyzer.effFile = cms.untracked.string( ivars.eff )
process.vnanalyzer.EPLevel = cms.untracked.int32(2)
process.vnanalyzer.Recenter = cms.untracked.bool(True)
process.vnanalyzer.minrun_ = cms.untracked.int32(326500)
process.vnanalyzer.maxrun_ = cms.untracked.int32(328500)
process.vnanalyzer.makeTree_ = cms.untracked.bool(False)
process.vnanalyzer.minvz_ = cms.untracked.double(-3.);
process.vnanalyzer.maxvz_ = cms.untracked.double(3.);
process.vnanalyzer.part_ = cms.untracked.string( ivars.part )
if ivars.part == 'LM':
    process.p = cms.Path(process.hltMB*process.offlinePrimaryVerticesRecovery*process.eventSelection*process.centralityBin*process.dbCent*process.QWEvent*process.QWV0Event*process.QWV0MVA*process.QWV0MVAVector*process.hiEvtPlane * process.hiEvtPlaneFlat*process.vnanalyzer)
if ivars.part == 'KS':
    process.p = cms.Path(process.hltMB*process.offlinePrimaryVerticesRecovery*process.eventSelection*process.centralityBin*process.dbCent*process.QWEvent*process.QWV0Event*process.QWV0MVA*process.QWV0MVAVector*process.hiEvtPlane * process.hiEvtPlaneFlat*process.vnanalyzer)
if ivars.part == 'CH':
    process.p = cms.Path(process.hltMB*process.offlinePrimaryVerticesRecovery*process.eventSelection*process.centralityBin*process.hiEvtPlane * process.hiEvtPlaneFlat*process.vnanalyzer)

from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
