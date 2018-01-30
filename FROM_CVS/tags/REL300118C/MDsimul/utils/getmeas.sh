#!/bin/sh
M=$HOME/simdat/measures
START_DIR=`pwd`
for dire in `ls -d -1 $1*` 
do
cd $M/$dire
echo -n -e $dire " "
#echo -n `echo $dire | awk '{ s=index("Vol",$0); e=index("NTV", $0); s=substr($0,s+1,e-s);\
# printf("%s",s(; }'`
$HOME/MDsimul/bin/md2ascii -t double -o Pex.a Pex-0
$HOME/MDsimul/utils/mean.sh Pex.a
rm Pex.a
done
cd $START_DIR
