if [ "$1" == "" ]
then
X0="1.5"
else
X0="$1"
fi
PERC="/home/demichel/HARDELL/Veff/X0_$X0"
NUMR="300"
IR="0"
ST="100000000"
TOTT="100000000000"
EXE="calc_veff_hardell"
while [ $IR -lt 300 ]
do
#calc_veff_hardell <trials> <X0> <NUMR> <ir> <savett>
cd ${PERC}/IR_$IR
EN="veff_IR_$IR"
ISRUN=`ps ax | grep $EN`
if [ "$ISRUN" == "" ]
then
echo "IR=" $IR " is not running, I will restart it"
nohup mosun ./$EN 0 >> screen_$IR &
fi
sleep 0.2
IR=$[$[IR]+1]
done
