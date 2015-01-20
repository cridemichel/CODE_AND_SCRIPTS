i=0
while [ $i -le 15 ]
do
echo "Doing " $i
echo `pwd`/ > lista.tmp
ls Cnf*_R$i >> lista.tmp
sort -t - -k 3 -n lista.tmp > lista_sorted.tmp
$HOME/MDsimul/bin/calcdrsq ./lista_sorted.tmp > drsq.dat_R$i
i=$[$i+1]
done
rm lista_sorted.tmp
rm lista.tmp
