#$1 = PN file $2=minimum $3 = file to average over
PNN="_histfixed.dat"
JF="_joined.dat"
if [ "$3" == "" ]
then
cat $1  | LC_NUMERIC=C gawk  -v min=$2  'BEGIN {s=0} {if ($1 < min){s+=$1*$2; norm+=$2;} else {s2+=$1*$2; norm2+=$2}} END {printf("iso dens=%f nem dens=%f\n",s/norm,s2/norm2}'
else
cat $1 | awk '{printf("%d %s\n",$1,$2)}' > $PNN
join -j 1 $PNN $3 > $JF
cat $JF | LC_NUMERIC=C gawk  -v min=$2  'BEGIN {s=0; s2=0} {if ($1 < min){s+=$3*$2; norm+=$2;} else {s2+=$3*$2; norm2+=$2}} END {printf("iso avg=%f nem avg=%f\n",s/norm,s2/norm2)}'
#rm -f $JF
#rm -f $PNN
fi
