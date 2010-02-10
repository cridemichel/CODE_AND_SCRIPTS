csize=0
#echo "NC=" $NC
ai=0  
#echo "COL1=" ${colarr[1]} "COL2=" ${colarr[2]} "NC=" $NC
#exit
head -1 $1 
cat $2 |
while read PRIMO RESTO
do
#echo "PRIMO=" $PRIMO "RESTO=" $RESTO> /dev/stderr
csizeold=$[csize]
csize=0
for par in $RESTO
do
csize=$[$csize+1]
done
#echo "CSIZE=" $csize  "CSIZEOLD=" $csizeold > /dev/stderr
if [ $[csize] -ne $[csizeold] ]
then
ai=$[$ai+1]
echo $ai >_ai_ 
#exit
fi
done
NUMCLS=`cat _ai_`
rm _ai_
echo "NUMERO CLUSTER=" $NUMCLS > /dev/stderr
#NUMCLS=$[$NUMCLS+1]
RED="1.0"
BLUE="0.0"
csize=0
FIRST=1
ai=0
cat $2 |
while read PRIMO RESTO
do
#echo " COL1=" $COL1 " COL2=" $COL2 " PRIMO=" $PRIMO " RESTO=" $RESTO
#echo $RESTO
csizeold=$[csize]
csize=0
#echo "REST=" $RESTO > /dev/stderr
#echo $RESTO | read -a oneline   
NP=0
for par in $RESTO
do
NP=$[$NP+1]
done
#HSL TO RGB
#H=`echo ${ai}/${NUMCLS}| bc -l`
if [ "$FIRST" == "1" ]
then
MAXSIZE=$NP
FIRST="0"
fi
#i cluster più grandi avranno hue=0 e saranno rossi
#mentre quelli più piccoli (NP=1) avranno hue dell'ordine
# di 0.7 che corrisponde al blue (basta vedere con gimp)
H=`echo "0.7*(1.0-$NP/$MAXSIZE)"| bc -l`
#H="0.75"
#echo "H=" $H "ai=" $ai "NUMCLS=" $NUMCLS > /dev/stderr
L="0.7"
S="0.8"
#S=`echo $ai/$NUMCLS |bc -l`
RED=`./hsl2rgb $H $S $L R 3`
GREEN=`./hsl2rgb $H $S $L G 3`
BLUE=`./hsl2rgb $H $S $L B 3`
#RED=`echo 1.0-$ai/$NUMCLS|bc -l`
#BLUE="0.0"
#GREEN=`echo $ai/$NUMCLS| bc -l`
for par in $RESTO
do
#continue
cat $1 | awk -v numpar=$par -v r=$RED -v g=$GREEN -v b=$BLUE 'BEGIN {fpA=0; fbB=0; np=0; cc=0}{if ($17=="C[red]") {if (np==numpar) fpA=1; np++}; if ($17=="C[green]") {if (np==numpar) fpB=1; np++} if (fpA==1) {if (cc==0) {for (i=1; i < 17; i++) printf("%s ",$i); printf("C[%s,%s,%s]\n",r,g,b); } else if (cc < 6) print $0; cc++; if (cc==6) {fpA=0;cc=0;}}; if (fpB==1) {if (cc==0) {for (i=1; i < 17; i++) printf("%s ",$i); printf("C[%s,%s,%s]\n",r,g,b);} else if (cc < 3) print $0; cc++; if (cc==3) {fpB=0;cc=0;} };}'
csize=$[$csize+1]
done
#echo "size=" $csize "sizeold=" $csizeold > /dev/stderr
if [ ! \( $[csize] -eq $[csizeold] \) ]
then
echo "NP="$NP "/"$MAXSIZE > /dev/stderr
echo "cluster index=" $ai  > /dev/stderr
ai=$[$ai+1]
fi
done
