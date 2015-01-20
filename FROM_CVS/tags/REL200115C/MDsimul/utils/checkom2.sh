#!/bin/sh
for f in T*/Om2*
do
A=`cat $f | awk -v fi=$f '{ if ($2 != "") print ("problems with ",fi) }'`
B=`wc -l $f | awk '{print $1}'`
if [ "$B" != "3000" ]
then
echo "wrong number of eigenvalues in " $f 
exit 1
fi
if [ "$A" != "" ]
then
echo "line with two values in " $f
exit 1
fi
done
