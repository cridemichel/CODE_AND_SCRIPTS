echo "Deleting unuseful Cnf files..."
if [ "$1" == "" ]
then
DIRS="T*"
else
DIRS=$1
fi
for f in $DIRS
do
cd  $f
echo "Dir: " $f
FN=`ls -1 Cnf* | sort -t - -k 3 -n | tail -1`
echo "file to save: " $FN
if [ ! -e tmp ]
then
mkdir tmp
fi
cp $FN tmp/
rm -f Cnf*
mv tmp/$FN .
cd ..
done
echo "done"
