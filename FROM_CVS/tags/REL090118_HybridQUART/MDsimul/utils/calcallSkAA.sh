i=0
while [ $i -le 15 ]
do
$HOME/MDsimul/bin/calcSkAA  Cnf*2000000*_R$i > SkAA_R$i.dat
i=$[$i+1]
done
