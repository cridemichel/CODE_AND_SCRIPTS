NUMR="300"
basename `pwd` > _aaa_
X0=`cat _aaa_ | awk -F _ '{print $2}'`
rm _aaa_
IR="0"
A="0.5"
B=`echo "${A}*${X0}" | bc -l`
DR=`echo 2.0*\(${B}-${A}\)/$NUMR | bc -l`
echo "A=" $A " B=" $B " DR=" $DR " X0=" $X0
while [ $IR -lt $NUMR ]
do
R=`echo "2.0*${A}+${DR}*${IR}+${DR}*0.5" |bc -l`
tail -1 IR_${IR}/veff_vs_tt.dat | awk -v r=$R '{printf("%.15G %.15G\n",r,$2)}'
IR=$[$[IR]+1]
done
