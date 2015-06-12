FN="super-story"
P1="part1"
P2="part2"
P3="part3"
echo -n > $P1
echo -n > $P2
echo -n > $P3
cc=0
FF=`head -n 1 lista`
echo "Primo file=" $FF
#NP=`cat $FF| awk '{if (NR==2) {print $0}; }'`
NP=`cat mols.dat | awk '{if (NR==1) {print $2};}'`
NPA=`cat mols.dat | awk '{if (NR==2) print $2}'`
BOX=`cat $FF| awk '{if (NR==6) {print $0}; }'`
echo "NP=" $NP " NPA=" $NPA
for f in `cat $1`
do
echo "file=" $f
cat $f | awk -v np=$NP -v npa=$NPA '{if (NR>319 && NR-319 <= npa) print ($3, $4, $5); }' >> $P2
cc=$[${cc}+1]
done
for f in `cat $1`
do
echo "file=" $f
cat $f | awk -v np=$NP -v npa=$NPA '{if (NR>319+npa && NR-319 <= np) print ($3, $4, $5); }' >> $P3
#cc=$[${cc}+1]
done
NPTOT=`echo "${cc}*${NP}" | bc`
NPATOT=`echo "${cc}*${NPA}" | bc`
NPSMALL=$[${NPTOT}-${NPATOT}]
echo "NPTOT= " $NPTOT " NPATOT= " $NPATOT " NPSMALL=" $NPSMALL
echo "Ntot $NPTOT" >  supermols.dat
#echo "Nlarge $NPATOT" >> supermols.dat
echo "0 0 " > $P1
echo $NPTOT >> $P1
echo $BOX >> $P1
echo "0 0 0" >> $P1
echo "0 0 0 0" >> $P1
cat $P1 $P2 $P3 > $FN
