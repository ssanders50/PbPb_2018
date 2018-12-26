import FWCore.ParameterSet.VarParsing as VarParsing
import string
import subprocess
import os

if os.access("flatparms_new.db",os.F_OK)==True:
    ret = subprocess.Popen(['rm','flatparms_new.db'])
    ret.wait()
    
ivars = VarParsing.VarParsing('standard')

ivars.register ('outputFile',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="for testing")

ivars.register ('outputTag',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="for testing")

ivars.register ('begin',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.int,
                info="for testing")

ivars.register ('end',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.int,
                info="for testing")


ivars.parseArguments()
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.CondDBCommon.connect = "sqlite_file:" + ivars.outputFile
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptyIOVSource",
    firstValue = cms.uint64(ivars.begin),
    lastValue = cms.uint64(ivars.end),
    timetype = cms.string('runnumber'),
    interval = cms.uint64(1)
)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
                                          process.CondDBCommon,
                                          timetype = cms.untracked.string("runnumber"),
                                          toPut = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                     tag = cms.string(ivars.outputTag)
                                                                     )
                                                            ),
                                          loadBlobStreamer = cms.untracked.bool(False)
                                          )


process.MoveFlatParamsToDB = cms.EDAnalyzer('MoveFlatParamsToDB',
                                           inputTFile = cms.string("data/rpflat_combined.root"),
                                           rootTag = cms.string(ivars.outputTag)
)


process.p = cms.Path(process.MoveFlatParamsToDB)
