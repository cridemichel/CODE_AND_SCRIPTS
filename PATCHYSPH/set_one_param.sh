#$1 = param file $2 = parametro $3 = valore 
cat $1 | awk -F ':' -v pn=$2 -v pv=$3 '{if ($1==pn) printf("%s: %s\n",pn, pv); else print $0;}' > PARFILE.TMP 
cp PARFILE.TMP $1
rm PARFILE.TMP
