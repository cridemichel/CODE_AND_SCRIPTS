PF="ellipsoid_flex_mc.par"
if [ ! -e tmp ]
then
mkdir tmp
fi
lastcnf=`ls -rt COORD_TMP_ASCII* CnfT* | tail -1`  
cd tmp
cp ../$lastcnf .
cat $lastcnf | awk 'BEGIN{at=0} {if (at >= 1) print $0; if ($1=="@@@") at++;}' > start.cnf
cp ../$PF .
cp ../../set_* ..
../set_params.sh $PF stepnum 1
EXE=`ls ../gapdna*`
$EXE -mgl 2 -fa $PF
