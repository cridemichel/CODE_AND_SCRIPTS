#!/bin/bash
cd /home/demichel/UNIVSCAL_POLYMERS/VISCO/
for j in `ls -d p_*`
do
ISR=`squeue --me | grep $j`
if [ "$ISR" == "" ]
then
echo "job $j is not running, I restarted it"
#batch /home/demichel/submit_$j
fi	
done
