#! /bin/bash
# aggpar.sh <regular expression per i files> <parametro> <nuovo valore> 
if [ "$1" == "" ]
then
echo "Nessun argomento, esco..."
exit
fi
for name in `ls $1`
do
echo "processing " $name
cp $name "TMP$name"
gawk -v par="$2" -v valo="$3" 'BEGIN {print (par ": " valo)}; \
  $1 == "@@@" {print ($0 "\n" (par ": " valo) }'\
 "TMP$name" > $name
rm "TMP$name"
done
