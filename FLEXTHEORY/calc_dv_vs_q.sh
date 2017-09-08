echo -n "" > v_vs_qsq.dat
echo -n "" > v_vs_q.dat
for f in `ls -d q_*`
do
QSQ=`echo $f| LANG=C awk -F _ '{print $2*$2}'`
Q=`echo $f| LANG=C awk -F _ '{print $2}'`
V=`tail -n 1 $f/v11.dat | LANG=C awk '{print $2}'` 
echo $QSQ $V >> v_vs_qsq.dat
echo $Q $V >> v_vs_q.dat
done
