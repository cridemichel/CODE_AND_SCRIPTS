#!/bin/sh
#$1 è un file contenente la lista delle configurazioni da distribuire sui nodi
i=0
DIR="$HOME/CRISTIANO/ANTIBODYGENE"
#bash array containing all available nodes
LN=("12" "15" "17" "7" "8" "26")
while [ $i -lt 6 ]
do 
NODE=${LN[$i]}
echo "NODE: irrmalin"${NODE}
ssh -n foffi@irrmalin${NODE} ${DIR}/kill_all_jobs_on_node.sh &
echo "JOBS KILLED" 
i=$[$i+1]
done
