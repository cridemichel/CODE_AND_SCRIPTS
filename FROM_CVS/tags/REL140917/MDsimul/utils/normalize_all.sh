#!/bin/sh
for $f in monohs_run?/*
do
NUM=`echo $f | awk -F 'run' '{print $2}' | awk -F "/" '{print $1}'` 
NOM=`echo $f | awk -F "/" {print $2}`
DST=${f}.norm$NUM
echo "processing" $f
cat $f | awk '{if (NR==1) scalfact=$2; print ($1,$2/scalfact)}' > $DST
done 

