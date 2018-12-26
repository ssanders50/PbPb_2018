from WMCore.Configuration import Configuration
config = Configuration()
from CRABClient.UserUtilities import getUsernameFromSiteDB
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.section_('JobType')
config.JobType.outputFiles = ['calib.root']
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/sanders/PbPb_2018/CalibTree/calibtree_cfg.py'
config.JobType.maxJobRuntimeMin = 1315
config.section_('Data')
config.Data.unitsPerJob = 80
config.Data.publication = False
config.Data.splitting = 'LumiBased'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_US_Vanderbilt'

config.Data.lumiMask = 'json_DCSONLY_HI.txt'
config.JobType.inputFiles = ['json_DCSONLY_HI.txt']

config.Data.runRange = '326381-326886'
config.General.requestName = 'PbPb2018_general_326381_326886'
config.Data.outLFNDirBase = '/store/user/ssanders/PbPb2018_general_326381_326886'
config.JobType.pyCfgParams = ['noprint','lumifile=json_DCSONLY_HI.txt','tracks=generalTracks']
config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v1/AOD'

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

submit(config)

config.Data.runRange = '326887-327564'
config.General.requestName = 'PbPb2018_general_326887_327564'
config.Data.outLFNDirBase = '/store/user/ssanders/PbPb2018_general_326887_327564'
config.JobType.pyCfgParams = ['noprint','lumifile=json_DCSONLY_HI.txt','tracks=generalTracks']
config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v2/AOD'
submit(config)

config.Data.runRange = '326381-326886'
config.General.requestName = 'PbPb2018_pixel_326381_326886'
config.Data.outLFNDirBase = '/store/user/ssanders/PbPb2018_pixel_326381_326886'
config.JobType.pyCfgParams = ['noprint','lumifile=json_DCSONLY_HI.txt','tracks=hiConformalPixelTracks']
config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v1/AOD'
submit(config)

config.Data.runRange = '326887-327564'
config.General.requestName = 'PbPb2018_pixel_326887_327564'
config.Data.outLFNDirBase = '/store/user/ssanders/PbPb2018_pixel_326887_327564'
config.JobType.pyCfgParams = ['noprint','lumifile=json_DCSONLY_HI.txt','tracks=hiConformalPixelTracks']
config.Data.inputDataset = '/HIMinimumBias0/HIRun2018A-PromptReco-v2/AOD'
submit(config)

