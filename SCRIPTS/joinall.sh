for f in EQUILB*
do 
echo "Doing " $f
cd $f 
if [ ! -e SHORT_TIMES ] 
then
cd ..
continue
fi
cd SHORT_TIMES
../../../joinsortSelf.sh
../../../joinsortColl.sh
../../../joinsortp.sh
cd ..
cd ..
echo "...done"
done
