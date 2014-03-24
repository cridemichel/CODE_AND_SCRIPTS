#!/bin/sh
if [ "$2" != "" ]
then
TIPO="$2"
else
TIPO="double"
fi
for ff in `ls $1` 
do
echo "processing " $ff "..."
$HOME/MDsimul/bin/md2ascii -t $TIPO -o $ff.dat $ff 
done
