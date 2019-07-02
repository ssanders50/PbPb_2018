import FWCore.ParameterSet.Config as cms

process = cms.Process("CumuV1")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.GeometryDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.HiGenCommon.AfterBurnerGenerator_cff')
process.load('PionExpoPt_cfi')



#process.AfterBurner = cms.Sequence(cms.SequencePlaceholder("AftBurner"))
process.AfterBurner = cms.Sequence(process.AftBurner)

process.RandomNumberGeneratorService.AftBurner = cms.PSet(
	initialSeed = cms.untracked.uint32(7410904),
	engineName = cms.untracked.string('HepJamesRandom')
)


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200000))
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options = cms.untracked.PSet(
		Rethrow = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source("EmptySource")

Mult = 500
part_id = cms.vint32();
for i in range(Mult):
        part_id.append(211)

process.generator.PGunParameters.PartID = part_id

#process.AftBurner.modv1 = cms.InputTag("0.")
process.AftBurner.modv1 = cms.InputTag("v1odd")
#process.AftBurner.modv2 = cms.InputTag("0.165646*exp(-( (x-2.64741)/1.36298 + exp( -(x-2.64741)/1.36298 ) )/2.)")
process.AftBurner.modv2 = cms.InputTag("0.0986502659")
process.AftBurner.fluct_v1 = cms.double(0.0)
process.AftBurner.fluct_v2 = cms.double(0.0)
process.AftBurner.modmethod = cms.int32(1)
#process.AftBurner.fixEP = cms.bool(False)



process.pgen = cms.Sequence(process.generator+cms.SequencePlaceholder("randomEngineStateProducer")+process.AfterBurner+process.GeneInfo)
#process.pgen = cms.Sequence(process.generator+cms.SequencePlaceholder("randomEngineStateProducer")+process.GeneInfo)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()


process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ampt_gen.root')
)

process.QWHIEvent = cms.EDProducer("QWHepMCProducer",
	src = cms.untracked.InputTag("generator")
)

process.QWGenEvent = cms.EDProducer("QWGenEventProducer",
	trackSrc  = cms.untracked.InputTag("genParticles"),
	ptMin = cms.untracked.double(0.3),
	ptMax = cms.untracked.double(3.),
	Etamin = cms.untracked.double(-5.0),
	Etamax = cms.untracked.double(5.0),
	isPrompt  = cms.untracked.bool(True)
)

process.AMPTCentrality = cms.EDProducer("QWIntProducer",
	src = cms.untracked.int32(61)
)

process.V1Cumu = cms.EDAnalyzer('QWV1Cumu'
	, trackEta = cms.untracked.InputTag('QWGenEvent', "eta")
	, trackPhi = cms.untracked.InputTag('QWGenEvent', "phi")
	, trackPt = cms.untracked.InputTag('QWGenEvent', "pt")
	, trackWeight = cms.untracked.InputTag('QWGenEvent', "weight")
	, trackCharge = cms.untracked.InputTag('QWGenEvent', "charge")
	, vertexZ = cms.untracked.InputTag('QWGenEvent', "vz")
	, RP = cms.untracked.InputTag('QWHIEvent', "EP")
	, centrality = cms.untracked.InputTag('AMPTCentrality')
	, minvz = cms.untracked.double(-15.0)
	, maxvz = cms.untracked.double(15.0)
	, rfpmineta = cms.untracked.double(-2.4)
	, rfpmaxeta = cms.untracked.double(2.4)
	, rfpminpt = cms.untracked.double(0.3)
	, rfpmaxpt = cms.untracked.double(3.0)
	, cmode = cms.untracked.int32(1)
	, nvtx = cms.untracked.int32(100)
)

process.histNoff = cms.EDAnalyzer('QWHistAnalyzer',
		src = cms.untracked.InputTag("AMPTCentrality"),
		Nbins = cms.untracked.int32(200),
		start = cms.untracked.double(0),
		end = cms.untracked.double(200),
		)

process.vectPhi = cms.EDAnalyzer('QWVectorAnalyzer',
		src = cms.untracked.InputTag("QWGenEvent", "phi"),
		hNbins = cms.untracked.int32(5000),
		hstart = cms.untracked.double(0),
		hend = cms.untracked.double(5000),
		cNbins = cms.untracked.int32(1000),
		cstart = cms.untracked.double(-3.14159265358979323846),
		cend = cms.untracked.double(3.14159265358979323846),
		)

process.vectPt = cms.EDAnalyzer('QWVectorAnalyzer',
		src = cms.untracked.InputTag("QWGenEvent", "pt"),
		hNbins = cms.untracked.int32(5000),
		hstart = cms.untracked.double(0),
		hend = cms.untracked.double(5000),
		cNbins = cms.untracked.int32(1000),
		cstart = cms.untracked.double(0),
		cend = cms.untracked.double(5),
		)

process.vectEta = cms.EDAnalyzer('QWVectorAnalyzer',
		src = cms.untracked.InputTag("QWGenEvent", "eta"),
		hNbins = cms.untracked.int32(5000),
		hstart = cms.untracked.double(0),
		hend = cms.untracked.double(5000),
		cNbins = cms.untracked.int32(48),
		cstart = cms.untracked.double(-5.),
		cend = cms.untracked.double(5.),
		)

process.RECO = cms.OutputModule("PoolOutputModule",
	outputCommands = cms.untracked.vstring('keep *'),
	SelectEvents = cms.untracked.PSet(
		SelectEvents = cms.vstring('ana')
		),
	fileName = cms.untracked.string('reco.root')
)

process.vectMon = cms.Sequence(process.histNoff*process.vectPhi*process.vectPt*process.vectEta)

process.ana = cms.Path(
		process.pgen *
		process.QWHIEvent *
		process.AMPTCentrality *
		process.QWGenEvent *
#		process.V1Cumu *
		process.vectMon
		)

#process.out = cms.EndPath(process.RECO)

#for path in process.paths:
#	getattr(process,path)._seq = process.generator * getattr(process,path)._seq
#
