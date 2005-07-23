PERC=$HOME/ELLIPSOIDS/FQT/
ls Cnf* > listaconf
q=2
while [ $q -lt 100 ]
do
$PERC/rho-dip << !
$q
!
q=$[$q+1]
done
