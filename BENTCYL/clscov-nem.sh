X0="2.0" # used only to overestimate box size here
OUTITS1="1000000"
OUTITSL="1000"
if [ \( "$1" == "" \) -o \( "$2" == "" \) -o \( "$3" == "" \) ]
then
echo "clscov.sh <seq> <size1> <size2> <alpha>"
echo "<seq> = ACC, AAC or AAT"
exit
fi
s1="$2"
s2="$3"
alpha="$4"
SEQN="$1"
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
L=`echo "2.0*${X0}*$pn+0.3*2.0*($s1+$s2-2)" | bc -l`
echo "$L $L $L" >> $fn
DN="alpha_$alpha"
if [ ! -e $DN ]
then
mkdir $DN
fi
cd $DN
pnm1=$[$s2]
if [ "$s2" == "1" ] 
then
OUTITS="$OUTITS1"
else
OUTITS="$OUTITSL"
fi
echo "100000000000 1 $pnm1 $OUTITS $alpha" > covmc.conf
#if [ ! -e $s1-$s2 ]
#then
#mkdir $s1-$s2
#fi
#cd $s1-$s2
#pnm1=$[$s2]
#echo "100000000000 0 $pnm1 5000" > covmc.conf
cp ../ellipsoid_flex.par .
../set_params.sh ellipsoid_flex.par parnum $pn inifile $fn
mv ../$fn .
EN="./bentHCCOV-$SEQN-$s1-$s2"
ln -sf ../bentHCCOV $EN
mosrun $EN -fa ellipsoid_flex.par > screen &
cd ..
