#!/bin/sh
awk 'BEGIN {tot=0; i=0} {tot += $2; i++} END {printf("\t\t#%d\t%.15G\n", i, tot/i);}' $1
