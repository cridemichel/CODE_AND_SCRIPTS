#!/bin/sh
if [ "$2" == "" ]
then 
P=100
else
P=$2
fi
for f in `ls $1`
do 
tail -n $P $f >> _mtail.tmp 
done
$HOME/MDsimul/utils/mean.sh _mtail.tmp
rm _mtail.tmp
