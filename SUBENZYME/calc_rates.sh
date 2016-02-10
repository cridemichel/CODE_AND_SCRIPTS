FN="k_vs_phi.dat"
FD="SP-popul.dat"
echo -n "" > $FN
for f in `ls -d -1 PHI_*`
do
cd $f
L=`tail -n 1 startEQ.cnf| awk '{print $1}'`
PN=`cat startEQ.cnf | grep parnum | awk -F : '{print $2}'`
PHI=`echo $f | awk -F _ '{print $2}'`
TIME=`tail -n 1 $FD | awk '{print $1}'`
EV=`tail -n 1 $FD | awk '{print $3}'`
echo "L= " $L " PN = " $PN
K=`echo "$EV/$TIME/($PN/$L/$L/$L)/2.0/2.0" | bc -l`
echo $PHI $K >> ../$FN
cd ..
done
