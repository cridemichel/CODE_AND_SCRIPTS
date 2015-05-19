EXE="CG-DNA-calc-v1-AF"
LISTA="7 9 12 20"
LISTAS="10 12 15 16"
INR="10"
for s in `echo $LISTAS`
do
if [ ! -e "s$s" ]
then
mkdir s$s
fi
cd s$s
for f in `echo $LISTA`
do
ir="1"
while [ $ir -le $INR ]
do
if [ ! -e "alpha_${f}_R$ir" ]
then
mkdir alpha_${f}_R$ir
fi
cd alpha_${f}_R$ir
ln -sf ../../$EXE
nohup mosrun ./$EXE ../../../CGDNA_AF/DD_${s}_PAS.pdb $s $f 10000000000 1 100000 > screen &
sleep 0.5
cd ..
#echo "$alpha $CV" >> $FN
ir=$[$ir+1]
done
done
cd ..
#s=$[$s+1]
done
