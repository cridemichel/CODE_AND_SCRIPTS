TEST=0
X0="0.4"
BASE=1.1
if [ "$1" == "" ]
then
echo "You must supply the volume fraction!"
exit
fi
i=1
while [ $i -le 10 ]
do
if [ ! -e X0_$X0 ]
then
mkdir X0_$X0
fi
echo "X0= "$X0 "Phi=" $1
if [ $TEST != "1" ]
then
../../SCRIPTS/sim1statepntPerf.sh $1 2 $X0
fi
i=$[$i+1]
X0=`echo "$X0*$BASE"|bc -l| awk '{printf("%1.3G",$1)}'`
done
