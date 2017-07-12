#!/bin/sh
# md2molgl.sh <regular expression per i files>
# TODO: this!!!! fare uno programma gawk a parte che scriva un file per molgl
# adeguato 
if [ "$1" == "" ]
then
echo "Nessun argomento, esco..."
exit
fi
for name in `ls $1`
do
echo "processing " $name
cp $name "$nameTMP"
gawk -f md2molgl.awk "$nameTMP" > $name
rm "$nameTMP"
done
