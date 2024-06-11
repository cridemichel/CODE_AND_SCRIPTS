#!/bin/bash
if [ "$1" = "" ]
then
	FR=0.8
else
	FR="$1"
fi
if [ "$2" == "" ]
then
	LD=`ls X0* -d`
else
	LD=$(cat $2)
fi
for f in $LD
do
	OD=`pwd`
	X0=$(echo $f | gawk -F _ '{print $2}')
	echo "Doing X0=" $f
	cd $f
	../ordparams.py $FR $X0 &
	cd $OD
done
