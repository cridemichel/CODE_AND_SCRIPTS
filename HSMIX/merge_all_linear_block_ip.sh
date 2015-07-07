P1="part1"
P2="part2"
P3="part3"
cc=0
BLOCK="$2"
FF=`head -n 1 lista`
echo "Primo file=" $FF
#NP=`cat $FF| awk '{if (NR==2) {print $0}; }'`
NP=`cat mols.dat | gawk '{if (NR==1) {print $2};}'`
NPA=`cat mols.dat | gawk '{if (NR==2) print $2}'`
BOX=`cat $FF| gawk '{if (NR==6) {print $0}; }'`
echo "NP=" $NP " NPA=" $NPA
#for f in `cat $1`
#do
#echo "file=" $f
#cat $f | awk -v np=$NP -v npa=$NPA '{if (NR>319 && NR-319 <= npa) print ($3, $4, $5); }' >> $P2
#cc=$[${cc}+1]
#done
aa=0
NFILES=`wc -l lista| gawk '{print $1}'`
while [ $cc -lt $NFILES ]
do
dd=0
FN="super-story-${aa}"
echo -n > $P1
echo -n > $P2
echo -n > $P3
num=`echo "$aa*$BLOCK+$BLOCK-1"|bc`
if [ $num -ge $NFILES ]
then
break
fi
while [ $dd -lt $BLOCK ]
do
num=`echo "$aa*$BLOCK+$dd"|bc`
f=`echo $num | gawk '{printf("story%04d",$1)}'`
echo "file=" $f
cat $f | gawk -v np=$NP -v npa=$NPA '{ip=NR-319; if (NR>319+npa && NR-319 <= np) print ($3, $4, $5, ip); }' >> $P3
dd=$[${dd}+1]
done
NPTOT=`echo "${BLOCK}*${NP}" | bc`
NPATOT=`echo "${BLOCK}*${NPA}" | bc`
NPSMALL=$[${NPTOT}-${NPATOT}]
echo "NPTOT= " $NPTOT " NPATOT= " $NPATOT " NPSMALL=" $NPSMALL
echo "Ntot $NPTOT" >  supermols.dat
#echo "Nlarge $NPATOT" >> supermols.dat
echo "$BLOCK 0" > $P1
echo $NPSMALL >> $P1
echo $BOX >> $P1
echo "0 0 0" >> $P1
echo "0 0 0 0" >> $P1
cat $P1 $P3 > $FN
cc=$[${cc}+${BLOCK}]
aa=$[${aa}+1]
done
