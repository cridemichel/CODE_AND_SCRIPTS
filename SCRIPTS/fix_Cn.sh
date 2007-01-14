for f in X0_*
do
cd $f
for ff in Phi*
do
cd $ff
if [ ! -e Cn.dat ]
then
cd ..
continue
fi 
NC=`cat Cn.dat | head -1 | awk '{print NF}'`
if [ "$NC" == "5" ]
then
echo -n "Processing " $f $ff " ..."
cp Cn.dat Cn1246.dat
cat Cn.dat | awk '{print ($1,$3,$4,$5)}' > Cn246.dat
mv Cn246.dat Cn.dat 
echo "done"
fi
cd ..
done
cd ..
done 

