FN="P_vs_phi.dat"
./calc_all_phi.sh
echo -n "" > $FN
for f in `ls -d P_*`
do
cd $f
NP=`wc -l phi.dat`
NEQ=`echo $NP/3| bc`
tail -n $NEQ phi.dat > _aaa_
AVPHI=`cat _aaa_ | awk '{c=c+1; sum+=$2} END {print sum/c}'`
rm _aaa_
PRESS=`echo $f|awk -F _ '{print $2}'`
echo $AVPHI $PRESS >> ../$FN
cd ..
done
