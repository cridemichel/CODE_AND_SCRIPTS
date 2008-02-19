declare -a colarr[0]
#echo "colori=" `cat $3`
read -a colarr < $3
NC=${#colarr[*]}
ai=0  
#echo "COL1=" ${colarr[1]} "COL2=" ${colarr[2]} "NC=" $NC
#exit
cat $2 |
while read PRIMO RESTO
do
if [ $[$ai+1] -ge $[NC] ]
then
COL1="red"
COL2="green"
else
COL1=${colarr[$ai]}
COL2=${colarr[$[$ai+1]]}
fi 
#echo " COL1=" $COL1 " COL2=" $COL2 " PRIMO=" $PRIMO " RESTO=" $RESTO
for par in $RESTO
do
cat $1 | awk -v numpar=$par -v col1=$COL1 -v col2=$COL2 'BEGIN {fpA=0; fbB=0; np=0; cc=0}{if ($17=="C[red]") {np++; if (np==numpar) fpA=1}; if ($17=="C[green]") {np++; if (np==numpar) fpB=1} if (fpA==1) {if (cc==0) {for (i=1; i < 17; i++) printf("%s ",$i); printf("C[%s]\n", col1); } else if (cc < 6) print $0; cc++; if (cc==6) {fpA=0;cc=0;}}; if (fpB==1) {if (cc==0) {for (i=1; i < 17; i++) printf("%s ",$i); printf("C[%s]\n",col2);} else if (cc < 3) print $0; cc++; if (cc==3) {fpB=0;cc=0;} };}'
done
ai=$[$ai+2]
done
