Name=$1
X0=$2

cat $Name |awk -v el=$X0 '{if(($2-el)==0.0) print $1, $3, $4}' |sort
