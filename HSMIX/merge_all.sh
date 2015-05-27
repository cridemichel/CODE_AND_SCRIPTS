FN="super-story"
P1="part1"
P2="part2"
P3="part3"
echo -n > $P2
cc=0
FF=`head -n 1 lista`
echo "Primo file=" $FF
NP=`cat $FF| awk '{if (NR==2) {print $0}; }'`
NPA=`cat mols.dat | awk '{if (NR==2) print $2}'`
BOX=`cat $FF| awk '{if (NR==3) {print $0}; }'`
echo "NP=" $NP " NPA=" $NPA
for f in `cat $1`
do
echo "file=" $f
cat $f | awk -v np=$NP -v npa=$NPA '{if (NR>5 && NR-5 <= npa) print $0; }' >> $P2
cc=$[${cc}+1]
done
for f in `cat $1`
do
echo "file=" $f
cat $f | awk -v np=$NP -v npa=$NPA '{if (NR>5+npa && NR-5 <= np) print $0; }' >> $P3
done
NPTOT=`echo "${cc}*${NP}" | bc`
NPATOT=`echo "${cc}*${NPA}" | bc`
echo "NPTOT= " $NPTOT " NPATOT= " $NPATOT
echo "Ntot $NPTOT" >  supermols.dat
echo "Nlarge $NPATOT" >> supermols.dat
echo "0 0 " > $P1
echo $NPTOT >> $P1
echo $BOX >> $P1
echo "0 0 0" >> $P1
echo "0 0 0 0" >> $P1
cat $P1 $P2 $P3 > $FN
