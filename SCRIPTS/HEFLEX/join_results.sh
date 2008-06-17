FN="bi-mono-bonds.dat"
echo -n "" > _lista_pref_
for f in `cat $1`
do 
PREF=`echo $f | awk -F R '{print $1}'`
echo $PREF >> _lista_pref_
done
for f in `cat _lista_pref_`
do
IGG=`echo $f | awk -F '/' '{print $1}'`
SIG=`echo $f | awk -F '/' '{print $2}'`
cat ${f}R*/$FN > joined_${IGG}_${SIG}_$FN
done
rm _lista_pref_
