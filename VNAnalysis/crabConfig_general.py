from WMCore.Configuration import Configuration
config = Configuration()
from CRABClient.UserUtilities import getUsernameFromSiteDB
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.section_('JobType')
config.JobType.outputFiles = ['vnanal.root']
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/sanders/PbPb_2018/VNAnalysis/vnanalysis_cfg.py'
config.JobType.maxJobRuntimeMin = 1350
config.section_('Data')
config.Data.allowNonValidInputDataset = True
config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v1/AOD'
config.Data.unitsPerJob = 50
config.Data.publication = False
config.Data.splitting = 'LumiBased'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_US_Vanderbilt'
#config.Site.storageSite = 'T2_US_MIT'
#config.Site.whitelist = ['T2_US_Vanderbilt']

config.Data.runRange = '326523-326886'
config.General.requestName = 'PbPb2018_vn_MB0'
config.Data.outLFNDirBase = '/store/user/ssanders/PbPb2018_vn_MB0'
config.Data.lumiMask = 'json_DCSONLY_HI.txt'
config.JobType.inputFiles=['json_DCSONLY_HI.txt','HeavyIonRPRcd_PbPb2018_general_offline.db','offset_PbPb2018_1_326886_general.root']
config.JobType.pyCfgParams = ['noprint','lumifile=json_DCSONLY_HI.txt','offset=offset_PbPb2018_1_326886_general.root','dbfile=HeavyIonRPRcd_PbPb2018_general_offline.db']
config.JobType.allowUndistributedCMSSW = True
if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
#print(" ================ ")
#print(config.Data.runRange)
#print(config.General.requestName)
#print(config.Data.outLFNDirBase)
#print(config.Data.lumiMask)
#print(config.JobType.inputFiles)
#print(config.JobType.pyCfgParams)
#submit(config)
iovs=[1, 326499, 326500, 326501, 326502, 326504, 326518, 326521, 326545, 326547, 326550, 326574, 326575, 326576, 326580, 326585, 
         326586, 326587, 326617, 326620, 326627, 326635, 326790, 326822, 326854, 326855, 326856, 326859, 326883, 326897, 326941, 
         326961, 326965, 326998, 327004, 327078, 327123, 327125, 327147, 327148, 327173, 327174, 327211, 327219, 327220, 327230, 
         327264, 327327, 327400, 327403, 327424, 328000]
niovs = len(iovs)-1;
for i in range(2, niovs):
    print(' =============== ')
    print(iovs[i])
    dataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v1/AOD'
    if iovs[i] > 326886:
        dataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v2/AOD'
    print(dataset)
    runrange = str(iovs[i])+'-'+str(iovs[i+1]-1)
    runranges = str(iovs[i])+'_'+str(iovs[i+1]-1)

    print(runrange)
    reqname = 'PbPb2018_vn_general_'+runranges
    print(reqname)
    dirbase='/store/user/ssanders/PbPb2018_vn_MB_general_'+runranges
    print(dirbase)
    infiles=['json_DCSONLY_HI.txt','HeavyIonRPRcd_PbPb2018_general_'+runranges+'.db','offset_PbPb2018_general_'+runranges+'.root']
    print(infiles)
    parms=['noprint','lumifile=json_DCSONLY_HI.txt','offset=offset_PbPb2018_general_'+runranges+'.root','dbfile=HeavyIonRPRcd_PbPb2018_general_'+runranges+'.db']
    print(parms)
    config.Data.inputDataset = dataset
    config.Data.runRange = runrange
    config.General.requestName = reqname
    config.Data.outLFNDirBase = dirbase
    config.JobType.inputFiles=infiles
    config.JobType.pyCfgParams = parms
    submit(config)
