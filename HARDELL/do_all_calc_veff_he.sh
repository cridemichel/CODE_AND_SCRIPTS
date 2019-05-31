NUMR="300"
IR="0"
X0="1.5"
ST="100000000"
TOTT="100000000000"
EXE="calc_veff_hardell"
while [ $IR -lt 300 ]
do
#calc_veff_hardell <trials> <X0> <NUMR> <ir> <savett>
mkdir IR_$IR
cd IR_$IR
cp ../$EXE .
EN="veff_IR_$IR"
ln -sf ./$EXE $EN
nohup mosrun ./$EN  $TOTT $X0 $NUMR $IR $ST > screen_$IR &
sleep 1.0
cd ..
IR=$[$[IR]+1]
done
