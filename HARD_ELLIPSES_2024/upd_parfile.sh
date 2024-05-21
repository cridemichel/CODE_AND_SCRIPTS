#!/bin/bash
# $1 = X0 directory 
# $2 = P  directory
FN="hellpars.xasc"
X0=$1
PS=$(echo $2| awk -F _ '{print $2}')
echo "X0=" $X0 " PS=" $PS
# reduced pressure Ps is \beta*P*2*a
# hence P = Ps/(2*a) since \beta=1
P=$PS # $(echo "$PS"| bc -l)
echo "[upd_parfile.sh] P=" $P " X0=" $X0
../../set_one_par.py $FN P $P 
../../set_one_par.py $FN a $X0
VAL=$(echo "2.0*$X0*1.01"| bc -l)
../../set_one_par.py $FN rcut $VAL
../../set_one_par.py $FN sigma $VAL
../../set_one_par.py $FN adjsteps 500000
