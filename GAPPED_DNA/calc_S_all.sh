for f in `ls -d P_*[0-9]`
do
cd $f
ls Cnf*| sort -t _ -k 3 -n > _lista_
NP=`wc -l _lista_ | awk '{print $1}'` 
NPR=`echo "$NP/3"| bc`
tail -n $NPR _lista_ > lista
rm _lista_
../order_param -mc -fl -c -t lista | awk '{print ($1, $4)}' > S.dat
cd ..
done
