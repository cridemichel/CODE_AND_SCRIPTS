#|/bin/sh
# scala <pattern dei file da scalare> <fattore> <tipo> <blocco>
# ved. anche scala.gwk per maggiori informazioni su tipo e blocco
if [ "$1" == "" ]
then
echo "Nessun argomento, esco..."
exit
fi
for name in `ls $1` 
do
gawk -v S=$2 -v tipo=$3 -v blocco=$4 -f ./scala.awk $name > "scaled-$name"
done
