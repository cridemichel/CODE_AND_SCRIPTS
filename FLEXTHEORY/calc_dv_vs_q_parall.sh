echo -n "" > v_vs_qsq.dat
echo -n "" > v_vs_q.dat
cd "alpha_20.0"
QS=`ls -d q_*`
cd ..
for f in `echo $QS` 
do
TTE="0.0"
TTR="0"
for a in `ls -d alpha_*`
do
cd $a
cd $f
TE=`tail -n 1 v11.dat | LANG=C awk '{print $3}'` 
TR=`tail -n 1 v11.dat | LANG=C awk '{print $1}'`
TTE=`echo ${TTE}+${TE} | bc -l`
TTR=`echo ${TTR}+${TR} | bc -l`
VOL3=`tail -n1 'iniconf_10-10.dat' |LANG=C awk '{print $1*$2*$3}'`
cd ..
cd ..
done
Q=`echo $f| LANG=C awk -F _ '{print $2}'`
QSQ=`echo $Q $Q| LANG=C awk '{print $1*$1}'`
EV=`echo ${TTE}*${VOL3}/$TTR| bc -l`
echo $QSQ $EV >> v_vs_qsq.dat
echo $Q $EV >> v_vs_q.dat
done
