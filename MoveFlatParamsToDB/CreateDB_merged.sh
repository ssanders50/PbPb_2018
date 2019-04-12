#!/bin/bash
#cd EPCalib
#rm -rf HiEvtPlaneList.h
#rm -rf HiEvtPlaneFlatten.h
#ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h
#ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h
#cd ..
rm -rf RescorTables
cp *.db dbfiles/.
rm macros/checkplots/merged*.pdf
rm tmpMerged.lis
ls -1 /panfs/crab_projects/crab_PbPb2018_merged_*/calib*.root > tmpMerged.lis
ls tmp.lis
declare -a list=(1 326545 326620 326887 327147 327230 328000)

nbreaks=${#list[@]}
echo $nbreaks
for (( i=0; i<${nbreaks}-1; i++));
do
min=${list[$i]}
max=${list[$i+1]}
max=$((max-1))
minrun=$min
maxrun=$max
range=$minrun
range+='_'
range+=$maxrun 
mkdir RescorTables
rm offsets/offset_PbPb2018_merged_$range.root
rm /rfs/sanders/tmp/EP_PbPb2018_merged_$range.root
rm  -rf RescorSave/RescorTables_PbPb2018_merged_$range_pixel
rm data/rpflat_combined.root
arg='EPCalib/EPCalib.C+('
arg+=$minrun
arg+=','
arg+=$maxrun
arg+=',"tmpMerged.lis","/panfs/tmpPbPb2018_merged","/panfs/EP_PbPb2018_merged_'
arg+=$range
arg+='.root","offsets/offset_PbPb2018_merged_'
arg+=$range
arg+='.root","RescorSave/merged_'
arg+=$range
arg+='")'
echo $arg
root -l -b -q $arg
ln -s /panfs/EP_PbPb2018_merged_$range.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_PbPb2018_merged_$range.db outputTag=HeavyIonRPRcd begin=$minrun end=$maxrun
rm /panfs/tmpPbPb2018_merged
mv RescorTables RescorSave/RescorTables_PbPb2018_merged_$range
conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_merged_$range.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_merged_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd -b $minrun -e $maxrun
cd macros
arg='CheckFlattening.C+("/panfs/EP_PbPb2018_merged_'
arg+=$range
arg+='.root")'
root -l -b -q $arg
cd ..
done
