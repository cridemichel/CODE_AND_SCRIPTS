PERC=$HOME/ELLIPSOIDS/FQT/
ls Cnf* | sort -t - -k 2 -k 3 -n > listaconf
q=2
while [ $q -lt 100 ]
do
$PERC/rho-dip << !
$q
!
q=$[$q+1]
done
