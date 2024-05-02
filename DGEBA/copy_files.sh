if [ "$1" = "" ]
then
echo "You must supply a Store file name"
exit 1
fi
FN="$1"
i=1
# CASO 2
while [ $i -le 20 ]
do
if [ ! -e RUN_$i ]
then
mkdir RUN_$i
fi
cp $HOME/UNIVSCAL_POLYMERS/CASO_2/RUN_${i}/$FN  RUN_${i}/store_file
cp $HOME/UNIVSCAL_POLYMERS/CASO_2/RUN_${i}/ellipsoid_dGEBA.par RUN_${i}/
i=$(($i+1))
done
# CASO 3
i=1
while [ $i -le 20 ]
do
j=$[$i+20]
if [ ! -e RUN_$j ]
then
mkdir RUN_$j
fi
cp $HOME/UNIVSCAL_POLYMERS/CASO_3/RUN_${i}/$FN  RUN_${j}/store_file
cp $HOME/UNIVSCAL_POLYMERS/CASO_3/RUN_${i}/ellipsoid_dGEBA.par RUN_${j}/
i=$(($i+1))
done
