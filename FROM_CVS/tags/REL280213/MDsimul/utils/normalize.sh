#!/bin/sh
awk '{ if (NR == 1) scalefact=$2; print ($1,$2/scalefact);}' $1 
