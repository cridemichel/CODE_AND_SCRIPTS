#!/bin/bash
DT="0.0001"
cc="0"
K="1.41421"
#TRUN="0.1"
TRUN="0.5"
OFL="logsigE_vs_logdt.dat"
OF="sigE_vs_dt.dat"
if [ \( "$1" == "triat" \) -o \( "$1" == "t" \) -o \( "$1" == "1" \) ]
then
PF="triatpars.xasc"
TPF="triatpars_templ.xasc"
EXE="../mdsim_triat 2"
else
PF="bmljpars.xasc"
TPF="bmljpars_templ.xasc"
EXE="../mdsim_bmlj 2"
fi
echo -n "" > $OF
echo -n "" > $OFL
while [ $cc -lt 14 ]
do
  cat $TPF >  $PF
  S=`echo "$TRUN/$DT" | bc -l | awk '{printf("%d",$0)}'`
  echo "dt: " $DT >> $PF
  echo "totsteps: " $S >> $PF
  cc=$[$cc+1]
  $EXE
  DE=`./calc_avg_sig.py energy.dat | awk '{print $2}'`
  AVGE=`./calc_avg_sig.py energy.dat | awk '{print $1}'`
  LDT=`echo "$DT" | awk '{print log($1)/log(10.0)}'`
  DEN=`echo "$DE $AVGE"| awk '{if ($2 < 0) ae=-$2; else ae=$2; printf("%.15G", $1/ae)}'`
  LDEN=`echo "$DEN" | awk '{printf ("%.15G", log($1)/log(10.0))}'`
  LDE=`echo $DE | awk '{printf("%.15G", log($1)/log(10.0))}'`
  echo $LDT $LDE >> $OFL
  echo $DT $DE >> $OF
  DT=`echo "$DT*$K" | bc -l`
done
