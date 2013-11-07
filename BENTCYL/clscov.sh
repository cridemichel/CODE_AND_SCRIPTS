X0="2.0" # used only to overestimate box size here
OUTITSL="10000"
OUTITS1="1000000"
if [ \( "$1" == "" \) -o \( "$2" == "" \) -o \( "$3" == "" \) ]
then
echo "clscov.sh <seq> <size1> <size2>"
echo "<seq> = ACC, AAC or AAT"
exit
fi
s1="$2"
s2="$3"
SEQN="$1"
pn=$[$s1+$s2]
fn="iniconf_${s1}-${s2}.dat"
cat conf_template_$SEQN > $fn
cat $fn | awk -v pn=$pn '{if ($0=="_PN1_") printf("parnum: %d\n",pn); else if ($0=="_PN2_") printf("%d\n",pn); else print $0}' > _aaa_
mv _aaa_ $fn
i=0
while [ $[i] -lt $[pn] ]
do
echo "0 0 0 1 0 0 0 1 0 0 0 1 0" >> $fn
i=$[$i+1]
done 
#big=`echo $s1 $s2| awk '{if ($1>$2) print $1; else print $2}'`
L=`echo "2.0*${X0}*$pn+0.15*2.0*($s1+$s2-2)" | bc -l`
echo "$L $L $L" >> $fn
if [ ! -e $s1-$s2 ]
then
mkdir $s1-$s2
fi
cd $s1-$s2
pnm1=$[$s2]
if [ "$s1" == "1" ]
then
echo "100000000000 0 $pnm1 $OUTITS1" > covmc.conf
else
echo "100000000000 0 $pnm1 $OUTITSL" > covmc.conf
fi
cp ../ellipsoid_flex_mc.par .
../set_params.sh ellipsoid_flex_mc.par parnum $pn inifile $fn
mv ../$fn .
EN="./bentHCCOV-$SEQN-$s1-$s2"
ln -sf ../bentHCCOV $EN
mosrun $EN -fa ellipsoid_flex_mc.par > screen &
cd ..
