i=0
while [ $i -le 15 ]
do
$HOME/MDsimul/bin/calcSkAAallk  Cnf*980000*_R$i > SkAAallk_R$i.dat
i=$[$i+1]
done
