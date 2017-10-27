#!/bin/bash
# getPE.sh <regexg per i files> <rango>
# <rango> si puo' non mettere se i files sono del tipo *_R<rango>
# dove <rango> è il rango del processo che ha scritto il file di restar processato
#echo -n  "" >  $2
if [ "$1" == "" ]
then
echo "Nessun argomento, esco..."
exit
fi
for name in `ls $1`
do
echo "processing " $name
gawk -v rango=$2 -v fn=$name -f ~/MDsimul/awksrc/getPE.awk "$name" > PE_of_$name.dat
#echo "&" > $2
done
