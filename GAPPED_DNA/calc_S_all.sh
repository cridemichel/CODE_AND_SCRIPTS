FNI2="I2_vs_phi.dat"
FNI4="I4_vs_phi.dat"
FNS="S_vs_phi.dat"
echo -n "" > $FNS
echo -n "" > $FNI2
echo -n "" > $FNI4
for f in `ls -d P_*[0-9]`
do
cd $f
ls Cnf*| sort -t _ -k 3 -n > _lista_
NP=`wc -l _lista_ | awk '{print $1}'` 
NPR=`echo "$NP/4"| bc`
tail -n $NPR _lista_ > lista
rm _lista_
../order_param_gapdna -mc -fl -c -t lista | awk '{print ($1, $4, $5, $6)}' > _aaa_
cat _aaa_ | awk '{print($1,$2)}' > S.dat
cat _aaa_ | awk '{print($1,$3)}' > I2.dat
cat _aaa_ | awk '{print($1,$4)}' > I4.dat
rm _aaa_
cd ..
done
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
