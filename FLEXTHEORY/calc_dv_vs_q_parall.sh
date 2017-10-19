echo -n "" > v_vs_qsq.dat
echo -n "" > v_vs_q.dat
cd "alpha_20.0_R0"
QS=`ls -d q_*`
TYPE="0"
if [ "$TYPE" == "0" ]
then
FN="covolume-nem.dat"
else
FN="v1${TYPE}.dat"
fi
cd ..
for f in `echo $QS` 
do
TTE="0.0"
TEX="0.0"
TTR="0"
for a in `ls -d alpha_*`
do
cd $a
cd $f
TE=`tail -n 1 $FN | LANG=C awk '{print $3}'` 
TEX=`tail -n 1 $FN | LANG=C awk '{print $4}'` 
TR=`tail -n 1 $FN | LANG=C awk '{print $1}'`
TTE=`echo ${TTE}+${TE} | bc -l`
TTEX=`echo ${TTEX}+${TEX}|bc -l`
TTR=`echo ${TTR}+${TR} | bc -l`
VOL3=`tail -n1 'iniconf_10-10.dat' |LANG=C awk '{print $1*$2*$3}'`
cd ..
cd ..
done
Q=`echo $f| LANG=C awk -F _ '{print $2}'`
QSQ=`echo $Q $Q| LANG=C awk '{print $1*$1}'`
EV=`echo (${TTE}-${TTEX})*${VOL3}/${TTR}| bc -l`
echo $QSQ $EV >> v_vs_qsq.dat
echo $Q $EV >> v_vs_q.dat
done
