cd $1
if [ ! -e tmp ]
then
mkdir tmp
fi
LC=`ls -rt COORD_TMP_ASCII*| tail -1`
EXE=`ls HC_N*| tail -1`
cd tmp
cp ../../ellipsoid_flex_mgl.par .
cp ../$LC .
cat $LC | awk 'BEGIN {at=0} {if (at>=1) print $0; if ($0=="@@@") at++;}' > startSUS.cnf
../$EXE -mgl 2 -fa ellipsoid_flex_mgl.par
cp startSUS.cnf.mgl ../conf.mgl
cd ..
#rm -fr tmp
cd ..
