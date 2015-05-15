ALFA="10 20 30"
for A in `echo $ALFA`
do 
FN="v2_vs_l_alpha$A.dat"
echo -n "" > $FN
DIRS="s1 s2 s3 s4 s5 s6 s7 s8 s9 s10"
for f in `echo $DIRS`
do
cd $f
cd alpha_$A
#echo "pwd=" `pwd`
AV=`tail -1 v2.dat| awk '{print $2/1.97^5/1E4}'`
S=`echo $f| awk -F s '{print $2}'`
echo "AV= " $AV " S=" $S
echo "$S $AV" >> ../../${FN}
cd ..  
cd ..
done
done
