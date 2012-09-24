if [ \( "$1" == "" \) -o \( "$2" == "" \) ]
then
echo "clscov.sh <size1> <size2>"
exit
fi
MOSRUN=""
OUTITS="50"
s1="$1"
s2="$2"
pn=$[$s1+$s2]
fn="iniconf_${s1}-${s2}.cnf"
cat conf_template > $fn
cat $fn | awk -v pn=$pn '{if ($0=="_PN1_") printf("parnum: %d\n",pn); else if ($0=="_PN2_") printf("%d\n",pn); else print $0}' > _aaa_
mv _aaa_ $fn
#PD is patch diameter
PD="0.119657"
#SIG is particle hard core diameter
SIG="1.0"
i=0
while [ $[i] -lt $[pn] ]
do
echo "0 0 0 1 0 0 0 1 0 0 0 1 0" >> $fn
i=$[$i+1]
done 
L=`echo "$SIG*$pn+($pn-2)*$PD*1.05"| bc -l`
echo "$L $L $L" >> $fn
if [ ! -e $s1-$s2 ]
then
mkdir $s1-$s2
fi
cd $s1-$s2
pnm1=$[$s2]
echo "100000000000 0 $pnm1 $OUTITS" > covmc.conf
cp ../ellipsoid_flex.par .
../set_params.sh ellipsoid_flex.par parnum $pn inifile $fn
mv ../$fn .
ln -sf ../ellipsCOV ellCOVlin-$s1-$s2
$MOSRUN ./ellCOVlin-$s1-$s2 -fa ellipsoid_flex.par > screen &
cd ..
