if [ "$1" = "" ]
then
echo "You must supply a par file name"
exit 1
fi
FN="$1"
i=1
while [ $i -le 20 ]
do
if [ ! -e RUN_$i ]
then
mkdir RUN_$i
fi
cp $HOME/UNIVSCAL_POLYMERS/CASO_2/RUN_${i}/$FN  RUN_${i}/
i=$(($i+1))
done
