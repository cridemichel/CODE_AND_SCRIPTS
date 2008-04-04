i=2
while [ $i -le 100 ]
do 
cat Fcoll-Q${i}.dat ../mediaColleCnfQ$i | sort -k 1 -n | awk '{printf("%.15G %.15G\n",0.75*$1,$2)}' > jointMediaColleCnfQ$i
i=$[$i+1]
done
