#!/bin/bash
for d in `ls -d RUN_*`
do
  OD=`pwd`
  cd $d
  ../../calcPACF.exe Ptens.dat 1000
  cd $OD
done	
