#!/bin/bash
for f in `ls X0_* -d`
do
	OD=`pwd`
	X0=$(echo $f | gawk -F _ '{print $2}')
	cd $f
	../eosth.py $X0
	cd $OD
done
