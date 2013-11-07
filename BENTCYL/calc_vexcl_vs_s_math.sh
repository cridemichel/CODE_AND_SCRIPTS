ALFA="5 10 15 20 25 30 35 40"
#ALFA="4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20"
FN="cov_vs_s_a.dat"
echo -n "{" > $FN
for A in `echo $ALFA`
do 
DIRS="s1 s2 s3 s4"
for f in `echo $DIRS`
do
cd $f
if [ -e alpha_$A ]
then
cd alpha_$A
echo "pwd=" `pwd`
AV=`tail -1 covolume-nem.dat| awk '{print $2}'`
S=`echo $f| awk -F s '{print $2}'`
echo "AV= " $AV " S=" $S
echo "{$S,$A,$AV}," >> ../../${FN}
cd ..  
fi
cd ..
done
done
echo "}" >> $FN
