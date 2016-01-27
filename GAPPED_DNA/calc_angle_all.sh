for f in P_*[0-9]
do
cd $f
cat dimers_info.dat| LANG=C gawk '{print $1, $3}' > avg-angle.dat
cd ..
done
FNA="angle_vs_phi.dat"
echo -n "" > $FNA
for f in `ls -d P_*[0-9]`
do
cd $f
echo "DIR=" $f
NP=`wc -l avg-angle.dat | awk '{print $1}'`
NEQ=`echo "${NP}"| bc`
echo "NP=" $NP " NEQ=" $NEQ
tail -n $NEQ avg-angle.dat > _aaa_
AVANG=`cat _aaa_ | awk '{c=c+1; sum+=$2} END {print sum/c}'`
rm _aaa_
NP=`wc -l phi.dat | awk '{print $1}'`
NEQ=`echo "${NP}"| bc`
echo "NP=" $NP " NEQ=" $NEQ
tail -n $NEQ phi.dat > _aaa_
AVPHI=`cat _aaa_ | awk '{c=c+1; sum+=$2} END {print sum/c}'`
rm _aaa_
PRESS=`echo $f|awk -F _ '{print $2}'`
echo $AVPHI $AVANG $PRESS >> ../$FNA
cd ..
done
