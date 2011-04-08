if [ \( "$1" == "" \) -o \( "$2" == "" \) ]
then
echo "clscov.sh <size1> <size2>"
exit
fi
s1="$1"
s2="$2"
pn=$[$s1+$s2]
fn="iniconf_${s1}-${s2}.dat"
cat conf_template > $fn
cat $fn | awk -v pn=$pn '{if ($0=="_PN1_") printf("parnum: %d\n",pn); else if ($0=="_PN2_") printf("%d\n",pn); else print $0}' > _aaa_
mv _aaa_ $fn
i=0
while [ $[i] -lt $[pn] ]
do
echo "0 0 0 1 0 0 0 1 0 0 0 1 0" >> $fn
i=$[$i+1]
done 
L=`echo 1.1*4.0*$pn | bc -l`
echo "$L $L $L" >> $fn
if [ ! -e $s1-$s2 ]
then
mkdir $s1-$s2
fi
cd $s1-$s2
pnm1=$[$s2]
echo "100000000000 0 $pnm1 5000" > covmc.conf
cp ../ellipsoid_flex.par .
../set_params.sh ellipsoid_flex.par parnum $pn inifile $fn
mv ../$fn .
ln -sf ../ellipsCOVMC ellCOVlin-$s1-$s2
mosrun ./ellCOVlin-$s1-$s2 -fa ellipsoid_flex.par > screen &
cd ..
