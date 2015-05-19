ALFA="7 9 12 20"
for A in `echo $ALFA`
do 
FN="v1_vs_l_alpha$A.dat"
echo -n "" > $FN
#DIRS="s1 s2 s3 s4 s5 s6 s7 s8 s9 s10"
DIRS="s10 s12 s15 s16"
for f in `echo $DIRS`
do
cd $f
#cd alpha_$A
#echo "pwd=" `pwd`
AV=`tail -q -n 1 alpha_${A}_R*/v1.dat| gawk 'BEGIN {tt=0; sum=0;} {tt+=$1; sum+=$2*$1} END {print sum/(1.97^4)/tt;}'`
S=`echo $f| awk -F s '{print $2}'`
echo "AV= " $AV " S=" $S
echo "$S $AV" >> ../${FN}
cd ..
done
done
