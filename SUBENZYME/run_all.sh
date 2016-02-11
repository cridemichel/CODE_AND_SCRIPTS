if [ "$1" == "" ]
then
LD=`ls -f PHI_*`
else
LD=`cat $1`
fi
SIG=`pwd | awk -F _ '{print $(NF-1)}'`
for f in `echo $LD`
do
cd $f
rm COORD_TMP[0,1]
cp CorFinal startEQ.cnf
PHI=`echo $f | awk -F _ '{print $2}'`
EN="subenz_sig_$SIG_phi_$PHI"
ln -sf ../subenz $EN
../set_params.sh ellipsoid_flex_prod.par stepnum 1000000 intervalSum 100.0
nohup mosrun ./$EN -fa ellipsoid_flex_prod.par > screen &
#mv CorFinal startGR.cnf
cd ..
sleep 0.5
done
