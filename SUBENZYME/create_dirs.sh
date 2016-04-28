if [ "$1" == "" ]
then
DIRS=`ls start_NOCROWD_MS*`
else
DIRS=`cat $1`
fi
for f in $DIRS
do
MS=`echo $f | awk -F _ '{print $4}'`
mkdir MS_${MS}
done
if [ "$1" == "" ]
then
#ln -sf ../grow_all.sh
ln -sf ../equilib_all.sh
ln -sf ../run_all.sh
cp ../set_* .
cp ../start.cnf .
cp ../ellipsoid_flex*par .
ln -sf ../subenz
fi
