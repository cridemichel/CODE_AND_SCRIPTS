Name=$1
Phi=$2

cat $Name |awk -v phi=$Phi '{if(($1-phi)==0.0) print $2, $3, $4}' |sort
