SC="$1"
NI="65278"
NE="338"
NC="160000"
DNS="3000"
DNS2="500"
SIGE="4.5"
SIGS="1.0"
SIGC="10.2"
Nav="6.02214129"
L=`tail -n 1 $1 | awk '{print $1}'`
NT=`echo "$NI+$NE+$NC"|bc`
i="0"
while [ $i -lt 22 ]
do
NS=`echo ${NI}-${i}*${DNS}|bc`
FF=`echo "1.0/0.1/${Nav}/${L}^3" | bc -l`
echo "NS=" $NS " FF=" $FF " L=" $L
MS=`echo "${FF}*${NS}"| bc -l | awk '{printf("%G",$1)}'`
echo "MS= " $MS
cat $SC | gawk -v nS=$NS -v nC=$NC -v nE=$NE -v nT=$NT  'BEGIN {at=0} {if (at==1 && NR==ATA+1) {NE=$1; NS=$2; NC=$4;}; if ($0=="@@@"){\
at=at+1; if (at==1) ATA=NR;\
else ATB=NR}; if (at==0 && $1=="parnum:") {print ("parnum:",nS+nE+nC);} else if (at==0) print $0; else \
if (at==1) {if (NR==ATA+1) print(nE, nS, "0", nC); else print $0;} else if (at==2) {if (NR-ATB<=nE) print $0; else if (NR-ATB<=nE+nS) print $0; else if (NR-ATB > NE+NS && NR-ATB <= nT) print $0; \
if (NR-ATB > nT) { if (NR-ATB <= nT+nE) print $0; else if (NR-ATB <= nT+nE+nS) print $0; else if (NR-ATB > nT+NE+NS && \
NR-ATB <= 2*nT) print $0;\
if (NR-ATB == 2*nT+1) print $0}}}' > start_CROWD_MS_${MS}_cnf
i=$[$i+1]
done
NI=$NS
#low concentrations with different DNS 
i="0"
while [ $i -lt 5 ]
do
NS=`echo ${NI}-${i}*${DNS2}|bc`
FF=`echo "1.0/0.1/${Nav}/${L}^3" | bc -l`
echo "NS=" $NS " FF=" $FF " L=" $L
MS=`echo "${FF}*${NS}"| bc -l | awk '{printf("%G",$1)}'`
echo "MS= " $MS
cat $SC | gawk -v nS=$NS -v nC=$NC -v nE=$NE -v nT=$NT  'BEGIN {at=0} {if (at==1 && NR==ATA) {NE=$1; NS=$2; NC=$4}; if ($0=="@@@"){\
at=at+1; if (at==1) ATA=NR;\
else ATB=NR}; if (at==0 && $1=="parnum:") {print ("parnum:",nS+nE);} else if (at==0) print $0; else \
if (at==1) {if (NR==ATA+1) print(nE, nS, "0", nC); else print $0;} else if (at==2) {if (NR-ATB<=nE) print $0; else if (NR-ATB<=nE+nS) print $0; else if (NR-ATB > NE+NS && NR-ATB <= nT) print $0; \
if (NR-ATB > nT) { if (NR-ATB <= nT+nE) print $0; else if (NR-ATB <= nT+nE+nS) print $0; else if (NR-ATB > nT+NE+NS && \
NR-ATB <= 2*nT) print $0;\
if (NR-ATB == 2*nT+1) print $0}}}' > start_CROWD_MS_${MS}_cnf
i=$[$i+1]
done
