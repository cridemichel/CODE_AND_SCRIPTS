i=0
while [ $i -le 15 ]
do
$HOME/MDsimul/bin/calcSkAB  Cnf*2000000*_R$i > SkAB_R$i.dat
i=$[$i+1]
done
