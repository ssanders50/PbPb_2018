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

ls -1 /rfs/sanders/crab_projects/crab_PbPb2018_MB0/calib*.root > tmp.lis
mkdir RescorTables
root -l  'EPCalib/EPCalib.C+(1,326886,"/rfs/sanders/tmp/tmpPbPb2018","/rfs/sanders/tmp/EPPbPb2018.root","/rfs/sanders/tmp/foffPbPb2018.root")'
rm offsets/offset_PbPb2018_1_326886.root
mv /rfs/sanders/tmp/foffPbPb2018.root offsets/offset_PbPb2018_1_326886_pixel.root
cd data
rm rpflat_combined.root
ln -s /rfs/sanders/tmp/EPPbPb2018.root rpflat_combined.root
cd ..
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_PbPb2018_1_326886.db outputTag=HeavyIonRPRcd begin=1 end=326886
rm /rfs/sanders/tmpPbPb2018
rm save/EP_PbPb2018_1_326886.root
mv /rfs/sanders/tmp/EPPbPb2018.root save/EP_PbPb2018_1_326886.root
rm  -rf RescorSave/RescorTables_PbPb2018_1_326886_pixel
mv RescorTables RescorSave/RescorTables_PbPb2018_1_326886_pixel
conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_1_326886.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_pixel_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_PbPb2018_offline -b 1 -e 326886

mkdir RescorTables
ls -1 /rfs/sanders/crab_projects/crab_PbPb2018_MB0_2/calib*.root > tmp.lis
root -l  'EPCalib/EPCalib.C+(326887,327489,"/rfs/sanders/tmp/tmpPbPb2018","/rfs/sanders/tmp/EPPbPb2018.root","/rfs/sanders/tmp/foffPbPb2018.root")'
rm offsets/offset_PbPb2018_326887_327489.root
mv /rfs/sanders/tmp/foffPbPb2018.root offsets/offset_PbPb2018_326887_327489_pixel.root
cd data
rm rpflat_combined.root
ln -s /rfs/sanders/tmp/EPPbPb2018.root rpflat_combined.root
cd ..
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_PbPb2018_326887_327489.db outputTag=HeavyIonRPRcd begin=326887 end=327489
rm /rfs/sanders/tmpPbPb2018
rm save/EP_PbPb2018_326887_327489.root
mv /rfs/sanders/tmp/EPPbPb2018.root save/EP_PbPb2018_326887_327489.root
rm  -rf RescorSave/RescorTables_PbPb2018_326887_327489
mv RescorTables RescorSave/RescorTables_PbPb2018_326887_327489

conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_326887_327489.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_pixel_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_PbPb2018_offline -b 326887 -e 327489
