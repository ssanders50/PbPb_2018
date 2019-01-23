#!/bin/bash
#cd EPCalib
#rm -rf HiEvtPlaneList.h
#rm -rf HiEvtPlaneFlatten.h
#ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h
#ln -s $CMSSW_BASE/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h
#cd ..
rm -rf RescorTables
cp *.db dbfiles/.
rm macros/checkplots/pixel*.pdf
rm tmpPixel.lis
ls -1 /rfs/sanders/crab_projects/crab_PbPb2018_pixel_326381_326886/calib*.root > tmpPixel.lis
ls -1 /rfs/sanders/crab_projects/crab_PbPb2018_pixel_326887_327564/calib*.root >> tmpPixel.lis
declare -a list=(1 326499 326500 326501 326502 326504 326518 326521 326545 326547 326550 326574 326575 326576 326580 326585 326586 326587 326617 326620 326627 326635 326790 326822 326854 326855 326856 326859 326883 326897 326941  326961 326965 326998 327004 327078 327123 327125 327147 327148 327173 327174 327211 327219 327220 327230 327264 327327 327400 327403 327424 328000)

nbreaks=${#list[@]}
echo $nbreaks
for (( i=42; i<${nbreaks}-2; i++));
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
rm offsets/offset_PbPb2018_pixel_$range.root
rm /rfs/sanders/tmp/EP_PbPb2018_pixel_$range.root
rm  -rf RescorSave/RescorTables_PbPb2018_pixel_$range_pixel
rm data/rpflat_combined.root
arg='EPCalib/EPCalib.C+('
arg+=$minrun
arg+=','
arg+=$maxrun
arg+=',"tmpPixel.lis","/rfs/sanders/tmp/tmpPbPb2018_pixel","/rfs/sanders/tmp/EP_PbPb2018_pixel_'
arg+=$range
arg+='.root","offsets/offset_PbPb2018_pixel_'
arg+=$range
arg+='.root","RescorSave/pixel_'
arg+=$range
arg+='")'
echo $arg
root -l $arg
ln -s /rfs/sanders/tmp/EP_PbPb2018_pixel_$range.root data/rpflat_combined.root
cmsRun moveflatparamstodb_cfg.py print outputFile=HeavyIonRPRcd_PbPb2018_pixel_$range.db outputTag=HeavyIonRPRcd begin=$minrun end=$maxrun
rm /rfs/sanders/tmpPbPb2018_pixel
mv RescorTables RescorSave/RescorTables_PbPb2018_pixel_$range
conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_pixel_$range.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_pixel_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_PbPb2018_offline -b $minrun -e $maxrun
cd macros
arg='CheckFlattening.C+("/rfs/sanders/tmp/EP_PbPb2018_pixel_'
arg+=$range
arg+='.root")'
root -l -b -q $arg
cd ..
done
