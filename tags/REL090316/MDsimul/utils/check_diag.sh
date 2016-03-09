#!/bin/sh
for om in `ls Om2*` 
do
head -n 3 $om | awk -v fi=$om 'BEGIN {ERR=0} {if ($1>0) v=$1; else v=-$1;\
 if (v >= 1E-8) ERR=1} END\
 {if (ERR) print ("autovalori negativi in " fi) } ' 
done
