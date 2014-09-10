L="101.234693"
XS=`"echo -$L/2"|bc -l`
YS=`"echo -$L/2"|bc -l`
DL="101.234693/16.0327"
NP=3050
for f in `cat $1`
do
i=1
j=1
xmax=$XS
ymax=$YS
while [ $i -lt 7 ]
do 
while [ $i -lt 7 ]
do
xmin=$xmax
ymin=$ymax
xmax=`echo $xmax+$DL| bc -l`
ymax=`echo $ymax+$DL| bc -l`
cat $f | awk -v np=$NP -v xm=$xmin -v ym=$ymin -v xM=$xmax -v yM=$ymax -v dl=$DL -f split_confs.awk 
i=$[$i+1]
done
j=$[$j+1]
done
done
