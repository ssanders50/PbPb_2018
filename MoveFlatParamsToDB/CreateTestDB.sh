#!/bin/bash
#cd EPCalib
#rm -rf HiEvtPlaneList.h
#rm -rf HiEvtPlaneFlatten.h
#ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h
#ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h
#cd ..
rm -rf RescorTables
cp *.db dbfiles/.
rm macros/checkplots/hi*.pdf
rm tmphi.lis
ls -1 /resfs/sanders/crab_projects/crab_PbPb2018_March2020_1_326544/calib*.root > tmphi.lis
ls tmphi.lis
declare -a list=(1 326545)

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
#mkdir RescorTables
rm offsets/offset_PbPb2018_March2020_$range.root
rm /resfs/sanders/tmp/EP_PbPb2018_March2020_$range.root
rm  -rf RescorSave/RescorTables_PbPb2018_March2020_$range_pixel
rm data/rpflat_combined.root
arg='EPCalib/EPCalib.C+('
arg+=$minrun
arg+=','
arg+=$maxrun
arg+=',"tmphi.lis","/resfs/sanders/tmpPbPb2018_March2020","/resfs/sanders/EP_PbPb2018_March2020_'
arg+=$range
arg+='.root","offsets/offset_PbPb2018_March2020_'
arg+=$range
arg+='.root","RescorSave/March2020_'
arg+=$range
arg+='")'
echo $arg
root -l -b -q $arg
ln -s /resfs/sanders/EP_PbPb2018_March2020_$range.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_PbPb2018_March2020_$range.db outputTag=HeavyIonRPRcd begin=$minrun end=$maxrun
rm /resfs/sanders/tmpPbPb2018_March2020
#mv RescorTables RescorSave/RescorTables_PbPb2018_March2020_$range
conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_March2020_$range.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_March2020_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd -b $minrun -e $maxrun
cd macros
arg='CheckFlattening.C+("/resfs/sanders/EP_PbPb2018_March2020_'
arg+=$range
arg+='.root")'
root -l -b -q $arg
cd ..
done
