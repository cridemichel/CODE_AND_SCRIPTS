#!/bin/bash
if [ "$1" = "" ]
then
	FR=0.5
else
	FR="$1"
fi
for f in `ls X0_* -d`
do
	OD=`pwd`
	X0=$(echo $f | gawk -F _ '{print $2}')
	cd $f
	../eos.py $FR $X0
	cd $OD
done
