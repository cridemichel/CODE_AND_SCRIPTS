#!/bin/sh
if [ "$1" == "" ]
then
PREPATH=""
else
PREPATH=$1
fi
for FF in `ls Qnc*_R*` 
do
if [ "$PREPATH" == "" ]
then
$HOME/MDsimul/bin/cri2fraMIX $FF ${FF}_FRA
else
$HOME/MDsimul/bin/cri2fraMIX $FF $PREPATH/${FF}_FRA
fi
done
