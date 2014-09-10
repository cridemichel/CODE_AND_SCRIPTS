L="101.234693"
XS=`echo -$L/2|bc -l`
YS=`echo -$L/2|bc -l`
DL="16.0327"
NP="3050"
cc=0
for f in `cat $1`
do
i=0
j=0
while [ $i -lt 5 ]
do 
while [ $j -lt 5 ]
do
xmin=`echo $XS+$DL*$i|bc -l`
ymin=`echo $YS+$DL*$j|bc -l`
xmax=`echo $XS+$DL*$i+$DL| bc -l`
ymax=`echo $YS+$DL*$j+$DL| bc -l`
cat $f | gawk -v Nigg=1000 -v Lbox=$L -v np=$NP -v xm=$xmin -v ym=$ymin -v xM=$xmax -v yM=$ymax -v dl=$DL -f split_confs.awk >_aaa_
RES=`tail -1 _aaa_`
if [ "$RES" != "FAILED" ]
then
NUMANT=`cat _NUMANT_`
NP=`echo $NUMANT+4|bc`
cat _aaa_ | awk -v numant=$NUMANT -v np=$NP '{if ($5=="_FIXNUMANT_") {print ("1 1 1 1", numant);} else if ($2=="_FIXPARNUM_") print("parnum:",np); else print $0; }' > monoGhostR$cc 
rm -f NUMANT 2>/dev/null
cc=$[$cc+1]
fi
#echo "BOH xmin= " $xmin, " xmax=", $xmax, " DL=" $DL
echo "i= ", $i, " j=" $j > /dev/stderr
j=$[$j+1]
done
j=0
i=$[$i+1]
done
done
