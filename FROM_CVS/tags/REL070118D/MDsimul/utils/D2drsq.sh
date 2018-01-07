#!/bin/sh
for fi in `ls D-0_R*.dat`
do
#~/MDsimul/bin/md2ascii -t double -o $fi ${fi}.dat
DRSQ=`echo "$fi" | awk -F _ '{print ("drSq-" $2)}'`  
echo $DRSQ
cat $fi | awk -v dt=$1 '{print ($1*dt "    "  $2*$1*dt*6)}' > $DRSQ
done
