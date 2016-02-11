for f in `ls -d sigE_*BIG`
do
cd $f
SIG=`echo $f | awk -F _ '{print $2}'`
ls -d -1 PHI_* > listadirs
echo "SIG= " $SIG
../expfit_scaled.sh listadirs $SIG
cd ..
done
