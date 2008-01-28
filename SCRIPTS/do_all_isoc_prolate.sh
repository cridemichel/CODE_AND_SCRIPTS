for f in Phi*
do
PHI=`echo $f | awk -F Phi '{print $2}'` 
cd $f
echo "Phi= " $PHI
../do_isocoreP.sh $PHI
cd ..
done
