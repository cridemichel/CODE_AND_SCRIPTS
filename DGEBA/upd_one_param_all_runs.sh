#!/bin/bash
if [ "$2" == "" ]
then
	echo "You must supply parametere name and new value"
	exit 1
fi
for d in `ls -d RUN_*`
do 
   OD=`pwd`
   cd $d
   ../../upd_param.sh $1 $2
   cd $OD
done
