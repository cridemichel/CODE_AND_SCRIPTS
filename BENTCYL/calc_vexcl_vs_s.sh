#ALFA="5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 25 30 35 40 45"
ALFA="5 10 15 20 25 30 35 40"
for A in `echo $ALFA`
do 
FN="cov_vs_s_a$A.dat"
echo -n "" > $FN
DIRS="s1 s2 s3 s4" #s5 s6"
for f in `echo $DIRS`
do
cd $f
cd alpha_$A
#echo "pwd=" `pwd`
AV=`tail -1 covolume-nem.dat| awk '{print $2}'`
S=`echo $f| awk -F s '{print $2}'`
echo "AV= " $AV " S=" $S
echo "$S $AV" >> ../../${FN}
cd ..  
cd ..
done
done
