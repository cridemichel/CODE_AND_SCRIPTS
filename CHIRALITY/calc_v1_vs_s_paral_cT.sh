ALFA="7 15 30"
for A in `echo $ALFA`
do 
FN="v1_vs_l_alpha$A.dat"
echo -n "" > $FN
#DIRS="s1 s2 s3 s4 s5 s6 s7 s8 s9 s10"
DIRS="s4 s6 s8 s10"
CONCS="500 600"
TEMPS="290 295"
for f in `echo $DIRS`
do
cd $f
#cd alpha_$A
#echo "pwd=" `pwd`
for c in `echo $CONS`
do
for t in `echo $TEMPS`
do
FNV="v1_c${c}_T${t}.dat"
AV=`tail -q -n 1 alpha_${A}_R*/$FNV| gawk 'BEGIN {tt=0; sum1=0; sum2=0; sum3=0} {tt+=$1; sum1+=$2*$1; sum2+=$3*$1; sum3+=$4*$1} END {print (sum1/tt, sum2/tt,sum3/tt);}'`
S=`echo $f| awk -F s '{print $2}'`
echo "AV= " $AV " S=" $S
echo "$S $AV" >> ../${FN}
done
done
cd ..
done
