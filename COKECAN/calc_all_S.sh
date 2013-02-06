echo -n  ""  > SvsN.dat
if [ "$1" == "" ]
then 
ls -d N_* | sort -t _ -k 2 -n > lista
LI="lista"
else
LI="$1"
fi
for f in `cat $LI`
do
cd $f
#ls Cnf* | sort -t _ -k 3 -n > lista
ls -rt COORD_TMP_ASCII* > lista
LST=`cat lista`
if [ "$LST" == "" ] 
then 
echo "no Cnf in dir " $f
cd ..
continue
fi
#ls -rt COORD_TMP_ASCII* > lista
if [ "$1" == "time" ]
then
../order_param -mc -fl -c -t lista > S.dat
else
NN=`echo $f| awk -F _ '{print $2}'` 
../order_param -mc -fl -c lista | awk -v nn=$NN '{printf("%d %f\n",nn,$3)}' >> ../SvsN.dat
fi
cd ..
done
