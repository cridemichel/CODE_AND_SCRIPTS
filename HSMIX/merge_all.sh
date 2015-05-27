FN="super-story"
P1="part1"
P2="part2"
echo -n > $P2
cc=0
FF=`head -n 1 lista`
echo "Primo file=" $FF
NP=`cat $FF| awk '{if (NR==2) {print $0}; }'`
BOX=`cat $FF| awk '{if (NR==3) {print $0}; }'`
for f in `cat $1`
do
echo "file=" $f
cat $f | awk '{if (NR==2) {NP=$0}; if (NR>5 && NR-5 <= NP) print $0; }' >> $P2
cc=$[${cc}+1]
done
echo "QUI cc=$cc NP=$NP"
NTOT=`echo $cc*$NP | bc`
echo "0 0 " > $P1
echo $NTOT >> $P1
echo $BOX >> $P1
echo "0 0 0" >> $P1
echo "0 0 0 0" >> $P1
cat $P1 $P2 > $FN
