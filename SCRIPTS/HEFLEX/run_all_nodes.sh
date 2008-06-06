#!/bin/sh
#$1 è un file contenente la lista delle configurazioni da distribuire sui nodi
i=0
LD=`cat $1`
DIR="$HOME/CRISTIANO/ANTIBODYGENE"
#bash array containing all available nodes
LN=("12" "12" "15" "15" "17" "17" "7" "7" "8" "8" "26" "26")
for CN in $LD
do 
NODE=${LN[$i]}
echo "NODE: irrmalin"${NODE}
echo "CONF: "${CN}  
ssh -n foffi@irrmalin${NODE} ${DIR}/run_one_job.sh ${CN} &
echo "JOB SUBMITTED" 
#ssh -n foffi@irrmalin${NODE} ${DIR}/renice_last_job.sh &
#echo "JOB RENICED"
i=$[$i+1]
done
