K1MAX="1.87"
K3MAX="2.17"
K1MAXB="3.97"
K3MAXB="3.22"
if [ "$2" != "" ]
then
K1MAX="$1"
K3MAX="$2"
fi
if [ "$4" != "" ]
then
K1MAXB="$3"
K3MAXB="$4"
fi
cat Welconst.dat | LANG=C awk -v K1M=$K1MAX -v K3M=$K3MAX '{if ($1 < K1M && $2 < K3M) print ($1,$2,$3)}' > W13.dat
cat Welconst.dat | LANG=C awk -v K1M=$K1MAXB -v K3M=$K3MAXB '{if ($1 < K1M && $2 < K3M) print ($1,$2,$4)}' > W23.dat
