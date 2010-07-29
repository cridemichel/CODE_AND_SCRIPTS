PF="silica_growth.par"
for f in PHI*
do
cd $f
PHI=`echo $f | awk -F _ '{print $2}'`
N=`echo $f | awk -F _ '{print $4}'`
COLL=`cat $PF | grep maxcoll | awk -F ':' '{print $2}'`
TI=`cat PRtime | grep real | awk '{print $2}'` 
#echo "COLL= " $COLL " TI" $TI
TC=`echo "${TI}/${COLL}"| bc -l`
echo "$N $PHI $TC" 
cd ..
done
