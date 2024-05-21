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
	cd $f
	../eos.py $FR
	cd $OD
done
