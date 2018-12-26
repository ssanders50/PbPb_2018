#!/bin/bash
cd EPCalib
rm *.so
rm *.d
rm *.pcm
rm -rf HiEvtPlaneList.h
rm -rf HiEvtPlaneFlatten.h
ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h
ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h
cd ..
rm -rf data/*.root
rm -rf RescorTables
rm *.db
rm /rfs/sanders/tmpPbPb2018
rm /rfs/sanders/EPPbPb2018.root
rm tmp.lis
rm macros/checkplots/general*.pdf

ls -1 /rfs/sanders/crab_projects/crab_PbPb2018_general_326381_326886/calib*.root > tmp.lis
mkdir RescorTables
rm offsets/offset_PbPb2018_general_1_326886.root
rm /rfs/sanders/tmp/EP_PbPb2018_general_1_326886.root
rm  -rf RescorSave/RescorTables_PbPb2018_general_1_326886_general
rm data/rpflat_combined.root
root -l  'EPCalib/EPCalib.C+(1,326886,"/rfs/sanders/tmp/tmpPbPb2018","/rfs/sanders/tmp/EP_PbPb2018_general_1_326886.root","offsets/offset_PbPb2018_general_1_326886.root")'
ln -s /rfs/sanders/tmp/EP_PbPb2018_general_1_326886.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_PbPb2018_general_1_326886.db outputTag=HeavyIonRPRcd begin=1 end=326886
rm /rfs/sanders/tmpPbPb2018
mv RescorTables RescorSave/RescorTables_PbPb2018_general_1_326886
conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_general_1_326886.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_general_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_PbPb2018_offline -b 1 -e 326886
cd macros
root -l -b -q 'CheckFlattening.C("/rfs/sanders/tmp/EP_PbPb2018_general_1_326886.root")'
cd ..

ls -1 /rfs/sanders/crab_projects/crab_PbPb2018_general_326887_327564/calib*.root > tmp.lis
mkdir RescorTables
rm offsets/offset_PbPb2018_general_326887_327564.root
rm /rfs/sanders/tmp/EP_PbPb2018_general_326887_327564.root
rm  -rf RescorSave/RescorTables_PbPb2018_general_326887_327564_general
rm data/rpflat_combined.root
root -l  'EPCalib/EPCalib.C+(1,327564,"/rfs/sanders/tmp/tmpPbPb2018","/rfs/sanders/tmp/EP_PbPb2018_general_326887_327564.root","offsets/offset_PbPb2018_general_326887_327564.root")'
ln -s /rfs/sanders/tmp/EP_PbPb2018_general_326887_327564.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_PbPb2018_general_326887_327564.db outputTag=HeavyIonRPRcd begin=326887 end=327564
rm /rfs/sanders/tmpPbPb2018
mv RescorTables RescorSave/RescorTables_PbPb2018_general_326887_327564
conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_general_326887_327564.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_general_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_PbPb2018_offline -b 326887 -e 327564
cd macros
root -l -b -q 'CheckFlattening.C("/rfs/sanders/tmp/EP_PbPb2018_general_326887_327564.root")'
cd ..

