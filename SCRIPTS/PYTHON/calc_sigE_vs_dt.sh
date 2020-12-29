DT="0.0001"
cc="0"
K="1.41421"
TRUN="0.5"
OF="sigE_vs_dt.dat"
PF="bmljpars.xasc"
TPF="bmljpars_templ.xasc"
echo -n "" > $OF
EXE="../mdsim_bmlj 2"
while [ $cc -lt 14 ]
do
  cat $TPF >  $PF
  S=`echo "$TRUN/$DT" | bc -l | awk '{printf("%d",$0)}'`
  echo "dt: " $DT >> $PF
  echo "totsteps: " $S >> $PF
  cc=$[$cc+1]
  $EXE
  DE=`./calc_avg_sig.py energy.dat | awk '{print $2}'`
  echo $DT $DE >> $OF
  DT=`echo "$DT*$K" | bc -l`
done
