FN="super-story"
P1="part1"
P2="part2"
echo -n > $P2
cc=0
NP=`cat `head -n 1 lista`| awk '{if (NR==2) {print $2}; }'`
BOX=`cat `head -n 1 lista`| awk '{if (NR==3) {print $0}; }'`
for f in `cat lista`
do
cat $f | awk '{if (NR==2) {NP=$2}; if (NR>5 && NR-5 <= NP) print $0 }' >> $P2
cc=$[${cc}+1]
done
NTOT=`echo $cc*$NP| bc`
echo "0 0 " > $P1
echo $NTOT >> $P1
echo $BOX >> $P1
echo "0 0 0" >> $P1
echo "0 0 0 0" >> $P1
cat $P1 $P2 > $FN
