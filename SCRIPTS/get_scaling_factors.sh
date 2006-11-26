echo -n "" > scaling_factors.dat 
for f in X0*
do
echo "Processing " $f
cd $f
PF=`ls Phi*/ellips.par -d | head -1` 
A0=`cat $PF | awk -F ':' '{if ($1=="A0") print $2}'`
B0=`cat $PF | awk -F ':' '{if ($1=="B0") print $2}'`
C0=`cat $PF | awk -F ':' '{if ($1=="C0") print $2}'`
#echo "ABC=" $A0 " " $B0 " " $C0
SF=`echo "e((1.0/3.0)*l($A0*$B0*$C0))"|bc -l`
AR=`echo $f| awk -F _ '{print $2}'`
cd ..
echo $AR " " $SF >> scaling_factors.dat
done
