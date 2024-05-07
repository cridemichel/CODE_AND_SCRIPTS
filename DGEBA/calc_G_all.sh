#!/bin/bash
for p in `ls -d p_*`
do
	PVAL=`echo $p | gawk -F _ '{print $2}'`
	tstar=$(cat pvststar.txt | gawk -v pval=$PVAL '{if (pval==$1){ print $2;}}')
	u2=$(cat pvststar.txt | gawk -v pval=$PVAL '{if (pval==$1){ print $3;}}')
        #echo "tstar=" $tstar
	OD=`pwd`
	cd $p
	G=$(cat media_acf_lin.dat | gawk -v ts=$tstar  'BEGIN {OLD=-1} {OLD=CUR; CUR=$1; if (ts >= OLD && ts < CUR) G=$2;  } END { print(G); }')	
	echo $u2 $G
	cd $OD
done
