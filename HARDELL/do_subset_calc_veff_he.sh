DIRS="IR_56 IR_122"
NUMR="300"
X0="1.1"
ST="1000000"
EXE="calc_veff_hardell"
for f in `echo $DIRS`
do
#calc_veff_hardell <trials> <X0> <NUMR> <ir> <savett>
cd $f
IR=`echo $f | awk -F _ '{print $2}'`
echo "IR=" $IR
nohup mosrun ./$EXE 100000000000 $X0 $NUMR $IR $ST > screen_$IR &
sleep 1.0
cd ..
done

