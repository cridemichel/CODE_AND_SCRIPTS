i=0
dr=0.5
FN="S_vs_raggio.dat"
echo "" > $FN
while [ $i -lt 35 ]
do
  ra=`echo "$i*$dr+4"| bc -l`
  echo "doing r=" $ra
  ./ordpar lista $ra | tail -n 1 >> $FN   
  i=$[$i+1]
done
