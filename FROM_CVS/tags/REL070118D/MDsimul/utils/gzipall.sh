if [ "$2" == "" ]
then
SUF='*FRA'
fi
for f in `ls -d $1` 
do
echo "gzipping" $f
gzip $f/$SUF
done
