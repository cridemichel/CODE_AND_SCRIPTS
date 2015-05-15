s="1"
EXE="CG-DNA-k2K22-v2-AF"
LISTA="10 20 30"
while [ $s -le 10 ]
do
if [ ! -e "s$s" ]
then
mkdir s$s
fi
cd s$s
for f in `echo $LISTA`
do
if [ ! -e "alpha_$f" ]
then
mkdir alpha_$f
fi
cd alpha_$f
ln -sf ../../$EXE
nohup mosrun ./$EXE ../../../CGDNA_AF/DD_${s}_PAS.pdb $s $f 10000000000 2 100000 > screen &
sleep 1
cd ..
#echo "$alpha $CV" >> $FN
done
cd ..
s=$[$s+1]
done
