i=1
while [ $[i] -le 4 ]
do
j=1
while [ $[j] -le 6 ]
do
if [ $[i] -gt $[j] ]
then
if [ ! -e ${i}-${j} ]
then
mkdir ${i}-${j}
fi
cp -f ${j}-${i}/covolume.dat ${i}-${j}/  
fi
j=$[$j+1]
done
i=$[$i+1]
done
