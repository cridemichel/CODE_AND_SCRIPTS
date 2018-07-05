for f in `cat $1`
do
NB=`echo $f | awk -F _ '{print $3}'| awk -F k '{print $2}'`
NQ=`echo $f | awk -F _ '{print $4}' | awk -F q '{print $2}'` 
FN="CNF-${NB}-${NQ}"
echo "NN: 38" > $FN
TI=`echo "1.3^40*$NB+1.3^$NQ"|bc -l`
echo "time: " $TI >> $FN
echo "@@@" >> $FN
echo "parnum: 56" >> $FN
echo "@@@" >> $FN
cat $f | awk '{print ($1,$2,"0.0")}'>> $FN
done
