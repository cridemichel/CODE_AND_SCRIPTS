P1="part1"
P2="part2"
P3="part3"
cc=0
#BLOCK="$2"
LST="$1"
FP="0"
LP="9"
FF=`head -n 1 $1`
echo "Primo file=" $FF
#NP=`cat $FF| awk '{if (NR==2) {print $0}; }'`
NP=`cat mols.dat | gawk '{if (NR==1) {print $2};}'`
NPA=`cat mols.dat | gawk '{if (NR==2) print $2}'`
BOX=`cat $FF| gawk 'BEGIN {r=-2} {if ($0=="A.SYSTEM SIZE") {r=NR;} if (NR==r+1) {print $0};}'`
L=`echo $BOX | awk '{print $1}'`
if [ $LP -gt $NP ]
then
LP="$NP"
fi
echo "NP=" $NP " NPA=" $NPA
#for f in `cat $1`
#do
#echo "file=" $f
#cat $f | awk -v np=$NP -v npa=$NPA '{if (NR>319 && NR-319 <= npa) print ($3, $4, $5); }' >> $P2
#cc=$[${cc}+1]
#done
aa=0
NFILES=`wc -l $1| gawk '{print $1}'`
dd=0
FN="super-story.mgl"
echo -n > $P1
echo -n > $P2
echo -n > $P3
num=`echo "$aa*$BLOCK+$BLOCK-1"|bc`
if [ $num -ge $NFILES ]
then
break
fi
#while [ $cc -lt $NFILES ]
for f in `cat $1`
do
# num=`echo "$aa*$BLOCK+$dd"|bc`
# f=`echo $num | gawk '{printf("story2-%04d",$1)}'`
#f=`echo $num | gawk '{printf("story%04d",$1)}'`
echo "file=" $f " FP" $FP " LP=" $LP
cat $f | gawk -v np=$NP -v lb=$L -v npa=$NPA -v fp=$FP -v lp=$LP 'BEGIN {REFNR=-1; co[0]="blue"; co[1]="green"; co[2]="cyan"; co[3]="orange"; co[4]="yellow"; co[5]="red"; co[6]="purple"; co[7]="turquoise"; co[8]="brown"; co[9]="green4"; } {if ($0=="H.LIST OF ATOMS") REFNR=NR+2; ip=NR-REFNR-npa-1; if (REFNR!=-1 && NR>=REFNR+npa+fp && NR-REFNR < lp+npa) {printf("%s %s %s @ %f C[%s]\n",$3-lb*0.5, $4-lb*0.5, $5-lb*0.5, 0.1, co[$1-npa-1]);} }' >> $P3
dd=$[${dd}+1]
done
NPTOT=`echo "${NFILES}*${NP}" | bc`
NPATOT=`echo "${NFILES}*${NPA}" | bc`
NPSMALL=$[${NPTOT}-${NPATOT}]
echo "NPTOT= " $NPTOT " NPATOT= " $NPATOT " NPSMALL=" $NPSMALL
echo "Ntot $NPTOT" >  supermols.dat
#echo "Nlarge $NPATOT" >> supermols.dat
VOL=`echo $L*$L*$L | bc -l`
echo ".Vol: " $VOL > $P1
#echo "$BLOCK 0" > $P1
#echo $NPSMALL >> $P1
#echo $BOX >> $P1
#echo "0 0 0" >> $P1
#echo "0 0 0 0" >> $P1
cat $P1 $P3 > $FN
#cc=$[${cc}+${BLOCK}]
#aa=$[${aa}+1]
