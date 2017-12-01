#!/bin/sh
$HOME/MDsimul/utils/bin2asc.sh 'D-0*'
$HOME/MDsimul/utils/D2drsq.sh $1
tar czf drsq.tgz drSq*dat
rm -f D-0*dat drSq*dat 
