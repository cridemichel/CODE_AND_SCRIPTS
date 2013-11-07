PR="./clscov.sh"
i=1
while [ $[i] -le 10 ]
do
echo "Launching $i-$i..."
$PR deg128 $i $i
sleep 1
i=$[$i+1]
done
