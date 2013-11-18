#!/bin/sh
# $1 = file Cnf
# $2 = numero particelle
if [ "$2" == "" ]
then
PAR=1000
else
PAR=$2
fi
if [ "$3" == "" ]
then
VOL=1000.0
else
VOL=$3
fi
echo ".Vol:" $VOL 
cat $1 | awk -v MOL=$PAR 'BEGIN { c=0; i=0; } { if (c == 2) { i+=1; if (i <=MOL) print $0;}; if ($1 == "@@@") c+=1;}'  
