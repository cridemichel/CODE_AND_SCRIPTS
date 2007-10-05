BP=$HOME/HDB/MD/MDSIMUL
ESE=$BP/bin/sticky
PARF="sticky.par"
CONV=$BP/fra2criSticky2-3
ST=SHORT_TIMES
if [ ! -e $ST ]
then
mkdir $ST 
fi
ls Cnf-*-0 | sort -t - -k 2 -n > lista
N=`wc -l lista | awk '{print $1}'`
i=0
echo "N=" $N
#while [ $i -lt $N ]
for f in `cat lista`
do
i=`echo $f | awk -F - '{print $2}'`
if [ ! -e $ST/CNF_$i ]
then
mkdir $ST/CNF_$i
fi
cp B*par $ST/CNF_$i/$PARF
T=`cat $ST/CNF_$i/$PARF | awk -F ':' '{if ($1=="temperat") print $2}'`
echo "T=" $T 
echo $T | $CONV Cnf-$i-0 $ST/CNF_$i/restart.res
#cat aaa | awk 'BEG {i=0} { if (i>0) print $0; if ($0=="@@@") i=i+1;}' > bbb
#cp aaa CNF_$i/restart.res
#rm -f aaa bbb
cd $ST/CNF_$i
../../../../set_one_param.sh $PARF storerate 0.01
../../../../set_one_param.sh $PARF stepnum 30
../../../../set_one_param.sh $PARF NN 30
$ESE -fa sticky.par
cd ..
cd ..
#i=$[$i+1]
done
