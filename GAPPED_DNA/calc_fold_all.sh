FN="foldfract_vs_phi.dat"
echo -n "" > $FN
./calc_phi_all.sh
for f in `ls -d P_*[0-9]`
do
cd $f
ls Cnf*| sort -t _ -k 3 -n > _lista_
NP=`wc -l phi.dat | awk '{print $1}'`
NEQ=`echo "${NP}/5"| bc`
echo "NP=" $NP " NEQ=" $NEQ
tail -n $NEQ phi.dat > _aaa_
AVPHI=`cat _aaa_ | awk '{c=c+1; sum+=$2} END {print sum/c}'`
NP=`wc -l _lista_ | awk '{print $1}'` 
NPR=`echo "$NP/5"| bc`
tail -n $NPR _lista_ > lista
rm _lista_
../calc_folding -mc -tr 30.0 lista 
FF=`cat foldfract.dat`
echo $AVPHI $FF >> ../$FN
cd ..
done
