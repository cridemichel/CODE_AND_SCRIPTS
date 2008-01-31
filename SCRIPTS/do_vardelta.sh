TEST=0
X0="$2"
BASE="1.165"
if [ "$1" == "" ]
then
echo "You must supply volume fractions and the elongation!"
exit
fi
for f in $1
do
VD="0.04"
cd Phi$f
echo "X0= "$X0 "Phi=" $f
i=1
while [ $i -le 20 ]
do
if [ $TEST != "1" ]
then
../sim1statepntPerfVardelta.sh $f 1 $X0 $VD
fi
i=$[$i+1]
VD=`echo "$VD*$BASE"|bc -l| awk '{printf("%1.5G",$1)}'`
done
cd ..
done
