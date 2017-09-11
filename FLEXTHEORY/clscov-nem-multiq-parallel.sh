X0="2.0"
QVECS="0.001 0.003 0.005 0.007 0.01 0.015 0.02 0.025 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1"
if [ \( "$1" == "" \) -o \( "$2" == "" \) ]
then
echo "clscov.sh <size1> <size2> <alpha>"
exit
fi
s1="$1"
s2="$2"
alpha="$3"
QVEC="$4"
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
NJOBS="9"
job="1"
while [ $[job] -le $[NJOBS] ]
do
echo "Running jobs #" ${job}
DN="alpha_${alpha}_R${job}"
if [ ! -e $DN ]
then
mkdir $DN
fi
cd $DN
for qv in `echo $QVECS`
do
mkdir "q_$qv"
cd q_$qv
pnm1=$[$s2]
echo "100000000000 11 $pnm1 5000 $alpha $qv" > covmc.conf
cp ../../ellipsoid_flex.par .
../../set_params.sh ellipsoid_flex.par parnum $pn inifile $fn
cp ../../$fn .
EXN="flexth-a${alpha}-$s1-$s2"
ln -sf ../../flextheory $EXN
mosrun ./$EXN -fa ellipsoid_flex.par > screen &
cd ..
done
cd ..
job=$[${job}+1]
done
