#!/bin/sh
rm eg.tmp
rm d2eg.tmp
if [ "$2" = "" ]
then
PTS=200
else
PTS=$2
fi
for f in $1
do
MEANEG=`$HOME/MDsimul/utils/meantails.sh $f $PTS| awk '{print $2}'`
echo "0 " $MEANEG >> eg.tmp
tail -n $PTS $f | awk '{print 0,$2}' | awk -v meg=$MEANEG 'BEGIN { co=0; sd2=0.0 } { co=co+1; sd2=sd2+($2-meg)*($2-meg)} END { print (0,sd2/co) }' >> d2eg.tmp
done
echo `$HOME/MDsimul/utils/mean.sh eg.tmp | awk '{print $2}'` " " `$HOME/MDsimul/utils/mean.sh d2eg.tmp | awk '{print $2}'`
