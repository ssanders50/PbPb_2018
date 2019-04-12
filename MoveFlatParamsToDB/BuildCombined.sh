
declare -a list=(1 326545 326620 326887 327147 327230 328000)
nbreaks=${#list[@]}
echo $nbreaks
for (( i=0; i<${nbreaks}-2; i++));
do
min=${list[$i]}
max=${list[$i+1]}
max=$((max-1))
minrun=$min
maxrun=$max
range=$minrun
range+='_'
range+=$maxrun 
conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_pixel_$range.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_pixel_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd_PbPb2018_offline -b $minrun -e $maxrun

conddb_import -f sqlite_file:HeavyIonRPRcd_PbPb2018_general_$range.db -c sqlite_file:HeavyIonRPRcd_PbPb2018_general_offline.db -i HeavyIonRPRcd -t HeavyIonRPRcd -b $minrun -e $maxrun
done
