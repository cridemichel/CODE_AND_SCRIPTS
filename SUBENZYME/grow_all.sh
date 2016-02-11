if [ "$1" == "" ]
then
LD=`ls -d PHI_*`
else
LD=`cat $1`
fi
SIG=`pwd| awk -F _ '{print $(NF-1)}'`
echo "SIG= " $SIG
for f in `echo $LD`
do
if [ ! -e $f ]
then
mkdir $f
fi
cd $f
rm COORD_TMP[0,1] 2> /dev/null
cp ../ellipsoid*par .
cp ../start.cnf .
PHI=`echo $f | awk -F _ '{print $2}'`
EN="subenz_sig_${SIG}_phi_$PHI"
ln -sf ../subenz $EN
../set_params.sh ellipsoid_flex_growth.par intervalSum 1.0 targetPhi $PHI
nohup mosrun ./$EN -fa ellipsoid_flex_growth.par > screen &
#mv CorFinal startGR.cnf
cd ..
sleep 0.5
done
