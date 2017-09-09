X0="2.0"
QCALC="0.08"
if [ \( "$1" == "" \) -o \( "$2" == "" \) ]
then
echo "clscov.sh <size1> <size2> <alpha>"
exit
fi
s1="$1"
s2="$2"
alpha="$3"
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
#big=`echo $s1 $s2| awk '{if ($1>$2) print $1; else print $2}'`
L=`echo "2.0*1.1*${X0}*$pn+0.16*2.0*($s1+$s2-2)" | bc -l | LANG=C gawk '{if ($1 < 4.0*1.1) {print 4.0*1.1} else {print $1}}'`
echo "$L $L $L" >> $fn
DN="alpha_$alpha"
if [ ! -e $DN ]
then
mkdir $DN
fi
cd $DN
mkdir "q_0"
cd "q_0"
pnm1=$[$s2]
echo "100000000000 1 $pnm1 5000 $alpha 0" > covmc.conf
cp ../../ellipsoid_flex.par .
../../set_params.sh ellipsoid_flex.par parnum $pn inifile $fn
mv ../../$fn .
EXN="flexth-a${alpha}-$s1-$s2"
ln -sf ../../flextheory $EXN
mosrun ./$EXN -fa ellipsoid_flex.par > screen &
cd ..
mkdir "q_$QCALC"
cd q_$CALC
pnm1=$[$s2]
echo "100000000000 11 $pnm1 5000 $alpha 0.08" > covmc.conf
cp ../../ellipsoid_flex.par .
../../set_params.sh ellipsoid_flex.par parnum $pn inifile $fn
mv ../../$fn .
EXN="flexth-a${alpha}-$s1-$s2"
ln -sf ../../flextheory $EXN
mosrun ./$EXN -fa ellipsoid_flex.par > screen &
cd ..
cd ..
