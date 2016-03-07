SC="$1"
NI="26662"
NE="338"
DNS="1500"
SIGE="4.5"
SIGS="1.0"
Nav="6.02214129"
L=`tail -n 1 $1 | awk '{print $1}'`
NT=`echo "$NI+$NE"|bc`
i="0"
while [ $i -lt 18 ]
do
NS=`echo ${NI}-${i}*${DNS}|bc`
FF=`echo "1.0/0.1/${Nav}/${L}^3" | bc -l`
echo "NS=" $NS " FF=" $FF " L=" $L
MS=`echo "${FF}*${NS}"| bc -l | awk '{printf("%G",$1)}'`
echo "MS= " $MS
cat $SC | gawk -v nS=$NS -v nE=$NE -v nT=$NT  'BEGIN {at=0} {if (at==1 && NR==ATA) {NE=$1; NS=$2}; if ("$0"=="@@@"){\
at=at+1; if (at==1) ATA=NR;\
else ATB=NR}; if (at==0 && "$1"=="parnum:") print ("parnum:",nS+nE); else if (at==0) print $0; else \
if (at==1) print $0; else if (at==2) {if (NR-(ATB+1)<nE) print $0; else if (NR-(ATB+1)<nE+nS) print $0; \
if (NR-(ATB+1) > nT) { if (NR-(ATB+1) < nT+nE) print $0; else if (NR-(ATB+1) < nT+nE+nS) print $0; \
if (NR-(ATB+1) == nT+1) print $0}}}' > start_NOCROWD_MS_${MS}.cnf
i=$[$i+1]
done
