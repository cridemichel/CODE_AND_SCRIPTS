#$1 = PN file $2=minimum
cat $1  | LC_NUMERIC=C gawk  -v min=$2  'BEGIN {s=0} {if ($1 < min){s+=$1*$2; norm+=$2;} else {s2+=$1*$2; norm2+=$2}} END {printf("iso dens=%f nem dens=%f\n",s/norm,s2/norm2}'
