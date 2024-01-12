#!/bin/bash
DT="0.0001"
cc="0"
K="1.41421"
#TRUN="0.1"
TRUN="1.0"
OFL="logsigE_vs_logdt.dat"
OF="sigE_vs_dt.dat"
if [ \( "$1" == "triat" \) -o \( "$1" == "t" \) -o \( "$1" == "1" \) ]
then
PF="triatpars.xasc"
TPF="triatpars_templ.xasc"
EXE="./mdsim_triat 2"
else
PF="bmljpars.xasc"
TPF="bmljpars_templ.xasc"
EXE="../mdsim_bmlj 2"
fi
echo -n "" > $OF
echo -n "" > $OFL
while [ $cc -lt 14 ]
do
  rm -f restart-*
  cat $TPF >  $PF
  S=`echo "$TRUN/$DT" | bc -l | awk '{printf("%d",$0)}'`
  #TAUP=`echo $DT 100.0| awk '{print ($1*$2)}'`
  #TAUB=`echo $DT 1000.0|awk '{print ($1*$2)}'`
  #echo "taup:" $TAUP >> $PF
  #echo "taub:" $TAUB >> $PF
  echo "dt: " $DT >> $PF
  echo "totsteps: " $S >> $PF
  cc=$[$cc+1]
  $EXE
  DE=`./calc_avg_sig.py energy.dat | awk '{print $2}'`
  AVGE=`./calc_avg_sig.py energy.dat | awk '{print $1}'`
  LDT=`echo "$DT" | LC=C gawk '{print log($1)/log(10.0)}'`
  DEN=`echo "$DE $AVGE"| LC=C gawk '{if ($2 < 0) ae=-$2; else ae=$2; printf("%.15G", $1/ae)}'`
  LDEN=`echo "$DEN" | LC=C gawk '{printf ("%.15G", log($1)/log(10.0))}'`
  LDE=`echo "$DE" | LC=C gawk '{printf("%.15G", log($1)/log(10.0))}'`
  echo $LDT $LDE >> $OFL
  echo $DT $DE >> $OF
  DT=`echo "$DT*$K" | bc -l | LC=C gawk '{printf("%f", $1)}'`
done
