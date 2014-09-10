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
while [ $i -lt 6 ]
do 
while [ $j -lt 6 ]
do
xmin=`echo $XS+$DL*$i|bc -l`
ymin=`echo $YS+$DL*$j|bc -l`
xmax=`echo $XS+$DL*$i+$DL| bc -l`
ymax=`echo $YS+$DL*$j+$DL| bc -l`
cat $f | LANG=C gawk -v Nigg=1000 -v Lbox=$L -v np=$NP -v xm=$xmin -v ym=$ymin -v xM=$xmax -v yM=$ymax -v dl=$DL -f split_confs.awk >_aaa_
RES=`tail -1 _aaa_`
if [ "$RES" != "FAILED" ]
then
NUMANT=`cat _NUMANT_`
NPR=`echo $NUMANT+4|bc`
cat _aaa_ | LANG=C gawk -v numant=$NUMANT -v np=$NPR '{if ($5=="_FIXNUMANT_") {print ("1 1 1 1", numant);} else if ($2=="_FIXPARNUM_") print("parnum:",np); else print $0; }' > _bbb_
#fix z coordinate
cat _bbb_ | LANG=C gawk -v Lbox=$L -v dl=$DL '{if (NF==13) print ($1,$2,$3+Lbox*0.5-dl*0.5,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13); else print $0;}' > _ccc_
# if any particles which belongs to IGG is out of the box along y pull it inside
DZ=`cat _ccc_ | LANG=C gawk -v dl=$DL 'BEGIN {dz=0} {if (NF==13 && $3 > dl*0.5 && $3-dl*0.5 > dz) dz=($3-dl*0.5); } END {print -dz*1.01}'`
cat _ccc_ | LANG=C gawk  -v dz=$DZ '{if (NF==13 && $13 <= 3) print ($1,$2,$3+dz,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13); else print $0}' > monoGhostSymFlex-0.200000-50.000000-$cc 
rm -f _NUMANT_ _aaa_ _bbb_ _ccc_ 2>/dev/null
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
