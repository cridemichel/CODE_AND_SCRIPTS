FN="super-story"
P1="part1"
P2="part2"
echo -n > $P2
cc=0
FF=`head -n 1 lista`
echo "Primo file=" $FF
NP=`cat $FF| awk '{if (NR==2) {print $0}; }'`
NPA=`cat mols.dat | awk '{if (NR==2) print $2}'`
BOX=`cat $FF| awk '{if (NR==3) {print $0}; }'`
echo "NP=" $NP " NPA=" $NPA
# only small particles
for f in `cat $1`
do
echo "file=" $f
cat $f | awk -v np=$NP -v npa=$NPA '{if (NR>5+npa && NR-5 <= np) print ("O ", $0);}' >> $P2
done
echo $NPTOT > $P1
echo  "superfile" ">> $P1
cat $P1 $P2 > $FN
rm -f $P1 $P2
