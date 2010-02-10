csize=0
declare -a colarr[0]
#echo "colori=" `cat $3`
read -a colarr < $3
NC=${#colarr[*]}
#echo "NC=" $NC
ai=0  
#echo "COL1=" ${colarr[1]} "COL2=" ${colarr[2]} "NC=" $NC
#exit
head -1 $1 
cat $2 |
while read PRIMO RESTO
do
if [ $[$ai+1] -ge $[NC] ]
then
echo "A) ai=" $ai "NC=" $NC > /dev/stderr
#default color if all grey levels are used
COL1="grey100"
COL2="grey100"
else
echo "B) ai=" $ai "NC=" $NC > /dev/stderr
COL1=${colarr[$ai]}
COL2=${colarr[$[$ai+1]]}
fi 
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
#NP=${#oneline[*]}
#monomeri, dimeri e trimeri casi speciali
echo "==> NP=" $NP > /dev/stderr
if [ $[NP] -eq 1 ]
then
COL1="green4"
COL2="green1"
elif [ $[NP] -eq 2 ]
then
COL1="blue4"
COL2="blue1"
elif [ $[NP] -eq 3 ]
then
COL1="magenta4"
COL2="magenta1"
fi
for par in $RESTO
do
#continue
cat $1 | awk -v numpar=$par -v col1=$COL1 -v col2=$COL2 'BEGIN {fpA=0; fbB=0; np=0; cc=0}{if ($17=="C[red]") {if (np==numpar) fpA=1; np++}; if ($17=="C[green]") {if (np==numpar) fpB=1; np++} if (fpA==1) {if (cc==0) {for (i=1; i < 17; i++) printf("%s ",$i); printf("C[%s]\n", col1); } else if (cc < 6) print $0; cc++; if (cc==6) {fpA=0;cc=0;}}; if (fpB==1) {if (cc==0) {for (i=1; i < 17; i++) printf("%s ",$i); printf("C[%s]\n",col2);} else if (cc < 3) print $0; cc++; if (cc==3) {fpB=0;cc=0;} };}'
csize=$[$csize+1]
done
echo "size=" $csize "sizeold=" $csizeold > /dev/stderr
if [ ! \( $[csize] -eq $[csizeold] \) ]
then
ai=$[$ai+2]
fi
done
