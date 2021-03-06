FNA="P_vs_phi.dat"
FNB="P_vs_S.dat"
./calc_phi_all.sh
echo -n "" > $FNA
echo -n "" > $FNB
for f in `ls -d P_*[0-9]`
do
cd $f
echo "DIR=" $f
NP=`wc -l phi.dat | awk '{print $1}'`
NEQ=`echo "${NP}/4"| bc`
echo "NP=" $NP " NEQ=" $NEQ
tail -n $NEQ phi.dat > _aaa_
AVPHI=`cat _aaa_ | awk '{c=c+1; sum+=$2} END {print sum/c}'`
NP=`wc -l S.dat | awk '{print $1}'`
NEQ=`echo "${NP}/4"| bc`
tail -n $NEQ S.dat > _aaa_
AVS=`cat _aaa_ | awk '{c=c+1; sum+=$2} END {print sum/c}'`
rm _aaa_
PRESS=`echo $f|awk -F _ '{print $2}'`
echo $AVPHI $PRESS >> ../$FNA
echo $AVS $PRESS >> ../$FNB
cd ..
done
