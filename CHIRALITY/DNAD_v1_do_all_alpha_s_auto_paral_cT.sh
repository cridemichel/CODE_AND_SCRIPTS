EXEN="CG-DNA-allAT-v1-P"
EXE="$HOME/CALCk2K22/allAT/v1_elec/PARALL/CG-DNA-allAT-v1-P"
LISTA="7 15 30"
LISTAS="4 6 8 10"
INR="10"
TEMP="temps.dat"
CONC="concs.dat"
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
ln -sf $EXE $EXEN
nohup mosrun $EXEN $HOME/CALCk2K22/allAT/CGDNA/allAT_${s}.pdb $s $f 10000000000 1 500000 500000 $TEMP $CONC 5.0 > screen &
sleep 1.1
cd ..
#echo "$alpha $CV" >> $FN
ir=$[$ir+1]
done
done
cd ..
#s=$[$s+1]
done
