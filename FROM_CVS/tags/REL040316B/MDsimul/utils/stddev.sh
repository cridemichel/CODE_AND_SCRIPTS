#!/bin/sh
rm _alleg_.tmp
for f in $1
do
$HOME/MDsimul/utils/meantails.sh $f 5 | awk '{print 0,$2}' >> _alleg_.tmp
done
MEANEG=`$HOME/MDsimul/utils/mean.sh _alleg_.tmp | awk '{print $2}'`
echo "media=" $MEANEG
cat _alleg_.tmp | awk -v meg=$MEANEG 'BEGIN { co=0; sd2=0.0 } { co=co+1; sd2=sd2+\
    ($2-meg)*($2-meg)} END { print (meg/9810.0,co,sd2/(co-1)/(9810.0*9810.0)) }' 
#rm _alleg_.tmp
