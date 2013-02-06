echo -n  ""  > ene_vs_N.dat
for f in N_*
do
cd $f
NN=`echo $f | awk -F _ '{print $2}'`
ENEM=`cat energy.dat | awk '{sum+=$2} END{print (sum/NR);}'`  
echo $NN $ENEM >> ../ene_vs_N.dat
cd ..
done
