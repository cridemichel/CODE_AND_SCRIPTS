X0="0.4"
BASE=1.1
if [ "$1" == ""]
then
echo "You must supply the volume fraction!"
exit
fi
i=1
while [ $i -le 10 ]
do
X0=`echo "$X0*$BASE"|bc -l`
if [ ! -e $X0 ]
then
mkdir $X0
fi
cd $X0
../../../SCRIPTS/sim1statepntPerf.sh $1 2 $X0
cd ..
i=$[$i+1]
done
