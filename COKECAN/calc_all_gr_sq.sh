TF="threshold.dat"
SQ="../calcSqDNAD"
GR="../calcgrDNAD"
LC="listacnf"
for f in PHI*
do
cd $f
if [ -e "$TF" ]
then
THR=`cat $TF`
else
THR="10000"
fi
ls Cnf*[^gz] | awk -F _ -v thr=$THR '{if ($3 > thr) print $0}' | sort -t _ -k 3 -n > $LC
$GR -nn -mc $LC 200
$SQ  -c $LC 0 150
cd ..
done
