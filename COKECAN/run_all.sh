PR="./clscov.sh"
i=1
while [ $[i] -le 4 ]
do
j=1
while [ $[j] -le 6 ]
do
if [ $[i] -eq $[j] ]
then
j=$[$j+1]
continue
fi
echo "Launching $i-$j..."
$PR $i $j
sleep 2
j=$[$j+1]
done
i=$[$i+1]
done
