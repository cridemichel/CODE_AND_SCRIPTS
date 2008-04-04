i=1
while [ $i -le 5 ]
do 
cat p${i}.dat ../p${i}.dat | sort -k 1 -n | awk '{printf("%.15G %.15G\n",0.75*$1,$2)}'> jointp$i.dat
i=$[$i+1]
done
