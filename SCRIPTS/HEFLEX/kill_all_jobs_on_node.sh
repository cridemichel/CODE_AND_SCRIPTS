#!/bin/sh
IDS=`ps ax | grep ellips-IgG | awk '{print $1}'` 
for f in $IDS
do
kill $f > /dev/null 2>&1
done
