i=0
rm lista
while [ $i -le 15 ]
do
echo `pwd`/ > lista.$i
ls Cnf*_R$i | sort -n -t - -k 4 >> lista.$i
$HOME/MDsimul/bin/buildfsk ./lista.$i > Fsk_R$i.dat
echo "Fsk_R$i.dat" >> filesselfCnf22.list 
i=$[$i+1]
done
$HOME/MDsimul/bin/mediafsk ./lista
