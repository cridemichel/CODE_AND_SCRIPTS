export LC_NUMERIC=C 
FN="bi-mono-bonds.dat"
if [ "$2" == "" ]
then
EQTIME="200000"
else
EQTIME="$2"
fi
if [ "$3" == "" ]
then
MAXTIME="1000000000"
else
MAXTIME="$3"
fi
if [ "$4" == "" ]
then
PERC="0.30"
else
PERC="$4"
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
if [ "$EQTIME" == "-1" ]
then
LV=`cat ${f}R*/$FN |tail -1| awk '{print $2}'` 
EQTIME2=`cat ${f}R*/$FN | awk -v lv="$LV" -v perc="$PERC" '{if ($2/lv-1 < perc && $2/lv-1 > -perc) {print $1; exit;}}'`
echo "LV=" $LV " EQTIME=" $EQTIME
else
EQTIME2="$EQTIME"
fi 
cat ${f}R*/$FN | awk -v eqtime="$EQTIME2" -v maxtime="$MAXTIME" '{if ($1 > eqtime && $1 < maxtime) print $0}'> joined_${IGG}_${SIG}_$FN
cd $IGG
if [ ! -e $SIGN ]
then
mkdir $SIGN
fi
cd $SIGN
mv ../../joined_${IGG}_${SIG}_$FN bi-mono-bonds.dat
i=0
while [ $i -lt 20 ]
do
if [ -e ../../${f}R${i} ]
then
cp ../../${f}R${i}/CorIni .
break
fi
i=$[$i+1]
done
echo $IGG/$SIGN/ >> ../../listaJoined
cd ..
cd ..
done
rm _lista_pref_
