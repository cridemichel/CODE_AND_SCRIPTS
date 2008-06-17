FN="bi-mono-bonds.dat"
if [ "$2" == "" ]
then
EQTIME="200000"
else
EQTIME="$2"
fi
echo -n "" > _lista_pref_
echo -n "" > listaJoined
for f in `cat $1`
do 
PREF=`echo $f | awk -F R '{print $1}'`
exist=0
for ff in `cat _lista_pref_`
do 
if [ "$ff" == "$PREF" ]
then
exist=1
fi
done
if [ "$exist" == "0" ]
then
echo $PREF >> _lista_pref_
fi
done
for f in `cat _lista_pref_`
do
IGG=`echo $f | awk -F '/' '{print $1}'`
SIG=`echo $f | awk -F '/' '{print $2}'`
SIGN=`echo $SIG | awk -F R '{print $1}'`
cat ${f}R*/$FN | awk -v eqtime="$EQTIME" '{if ($1 > eqtime) print $0}'> joined_${IGG}_${SIG}_$FN
cd $IGG
if [ ! -e $SIGN ]
then
mkdir $SIGN
fi
cd $SIGN
mv ../../joined_${IGG}_${SIG}_$FN bi-mono-bonds.dat
cp ../../${f}R1/CorIni .
echo $IGG/$SIGN/ >> ../../listaJoined
cd ..
cd ..
done
rm _lista_pref_
