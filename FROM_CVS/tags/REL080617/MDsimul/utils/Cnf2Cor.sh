#!/bin/sh
# $1 = file Cnf
# $2 = numero particelle
cat $1 | awk  'BEGIN { c=0; } { if (c > 0) {print $0;}; if ($1 == "@@@") c=1;}'  
