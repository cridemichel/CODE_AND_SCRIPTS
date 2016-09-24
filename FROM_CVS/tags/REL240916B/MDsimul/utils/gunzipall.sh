if [ "$2" == "" ]
then
SUF='*gz'
fi
for f in `ls -d $1` 
do
echo "gUNzipping" $f
gzip -d $f/$SUF
done
