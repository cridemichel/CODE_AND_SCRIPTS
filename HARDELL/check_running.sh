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
STATUS=`cat ${PERC}/CHECK_STATUS`
if [ \( "$STATUS" == "OFF" \) -o \( "$STATUS" == "off" \) -o \( "$STATUS" == "0" \) ]
then
echo "Check disabled..."
exit
fi
while [ $IR -lt 300 ]
do
#calc_veff_hardell <trials> <X0> <NUMR> <ir> <savett>
if [ ! -e ${PERC}/IR_$IR ]
then
continue
fi
cd ${PERC}/IR_$IR
EN="veff_IR_$IR"
ISRUN=`ps ax | grep $EN| grep mosrun`
if [ "$ISRUN" == "" ]
then
echo "IR=" $IR " is not running, I will restart it"
if [ -e "calcveff.chk" ]
then
nohup mosrun ./$EN $TOTT $X0 $NUMR $IR $ST > screen_$IR &
else 
nohup mosrun ./$EN 0 >> screen_$IR &
fi
fi
sleep 0.2
IR=$[$[IR]+1]
done
