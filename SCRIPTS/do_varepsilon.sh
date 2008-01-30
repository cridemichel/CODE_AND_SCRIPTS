TEST=1
X0="$2"
VE="0.0000001"
VEF="0.00000011"
BASE="1.4126"
if [ "$1" == "" ]
then
echo "You must supply volume fractions and the elongation!"
exit
fi
for f in $1
do
cd Phi$f
echo "X0= "$X0 "Phi=" $f
i=1
while [ $i -le 20 ]
do
echo "epsilon= " $VE "epsilonF= " $VEF
if [ $TEST != "1" ]
then
../sim1statepntPerfVarepsilon.sh $f 1 $X0 $VE $VEF
fi
i=$[$i+1]
VE=`echo "$VE*$BASE"|bc -l| awk '{printf("%1.5G",$1)}'`
VEF=`echo "$VE*$BASE"|bc -l| awk '{printf("%1.5G",$1)}'`
done
cd ..
done
