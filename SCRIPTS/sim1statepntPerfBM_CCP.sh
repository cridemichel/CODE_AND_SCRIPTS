# $1 = number of particles
# $2 = total steps 
# $3 = equilibration steps
# $4 = Q (size ratio big/small)
# $5 = extra label
# $6 = USEMLL
# $7 = pre-growth steps (if 0 then do not grow)
# $8 = volume fraction 
# 
if [ "$1" = "" ]
then
echo "Syntax: sim1statepnt <initial conf> <steps> <equilibration steps> <rcut>"
exit 0
fi
if [ "$2" = "" ]
then
echo "You have to supply total steps"
exit 0
fi
TOTSTP="$2"
if [ "$3" = "" ]
then
echo "You have to supply equilibration steps"
exit 0
fi
EQSTP="$3"
if [ "$4" = "" ]
then
Q="1.0"
else
Q="$4"
fi
if [ "$5" = "" ]
then
EXTLAB="BMHS"
else
EXTLAB="$5"
fi
GENF="./genconfHS"
#we fix number of big particles
TEMP="1.0"
TIMEEXE="/usr/bin/time"
PARFILE="silica_growth.par"
INTSUM="0.2"
INTSUMGR="0.05"
if [ "$6" == "" ]
then
USEMLL="1"
else
USEMLL="$0"
fi
if [ "$7" == "" ]
then
GRSTP="0"
else
GRSTP="$7"
fi
if [ "$8" == "" ]
then
PHI=`echo $INIFILE | awk -F 'Phi' '{print $2}'| awk -F _ '{print $1}'`
else
PHI="$8"
fi
SETPARAMS="../../set_params.sh"
PD="PHI_${PHI}"
if [ ! -e $PD ]
then
mkdir $PD
fi
N1=250
N2=`echo "$1-$N1" | bc`
INIFILE="start.cnf"
$GENF $N1 $N2 > $INIFILE
cp $PARFILE $PD
cp $INIFILE $PD
cd $PD
rm -f COORD_TMP*
rm -f Store-*
if [ "$USEMLL" == "1" ]
then 
BMHSEXE="../stickyMLL"
else
BMHSEXE="../stickyLL"
fi
SIMPR="bmhs_N_$1_Q_$Q_Phi${PHI}_${EXTLAB}_PR"
SIMEQ="bmhs_N_$1_Q_$Q_Phi${PHI}_${EXTLAB}_EQ"
STORERATE="0.0"
#PARNUM=512
#PARNUMA=512
DT="0.05"
RNNL="0.12"
#
#
#growth run
if [ "$GRSTP" != "0" ]
then
$SETPARAMS $PARFILE useNNL $USENNL Dt $DT stepnum $GRSTP VSteps 0 temperat $TEMP scalevel 1 rescaleTime 0.5 targetPhi 0.0 storerate $STORERATE intervalSum $INTSUM DtrCalc 0 rcut $RCUT rotMSDCalc 0 rmsd2end -1.0 tmsd2end -1.0 inifile $INIFLOC endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMEQ
./$SIMEQ -fa $PARFILE > screen_$SIMEQ 
#
$SETPARAMS $PARFILE useNNL $USENNL Dt $DT stepnum 500000000 VSteps 0 temperat $TEMP scalevel 0 rescaleTime 0.0 targetPhi $PHI storerate $STORERATE intervalSum $INTSUMGR DtrCalc 0 rcut $RCUT rotMSDCalc 0 rmsd2end -1.0 tmsd2end -1.0 inifile ${SIMEQ}.cor endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMEQ
./$SIMEQ -f $PARFILE > screen_$SIMEQ 
INIFEQ="CorFinal"
else
INIFEQ=$INIFLOC
fi
#else
#INIFEQ="$INIFLOC"
#fi
#equilibration run
$SETPARAMS $PARFILE useNNL $USENNL Dt $DT stepnum $EQSTP VSteps 0 temperat $TEMP scalevel 1 rescaleTime 0.5 targetPhi 0.0 storerate $STORERATE intervalSum $INTSUM DtrCalc 0 rcut $RCUT rotMSDCalc 0 rmsd2end -1.0 tmsd2end -1.0 inifile $INIFEQ endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMEQ
./$SIMEQ -fa $PARFILE > screen_$SIMEQ 
#
#
#production run
$SETPARAMS $PARFILE useNNL $USENNL stepnum $TOTSTP scalevel 0 Steps 0 temperat $TEMP targetPhi 0.0 storerate $STORERATE intervalSum $INTSUM DtrCalc 0 rcut $RCUT rotMSDCalc 0 rmsd2end -1.0 tmsd2end -1.0 inifile ${SIMEQ}.cor endfile ${SIMPR}.cor
ln -sf $ELLEXE $SIMPR
if [ -e /Applications ]
then
$TIMEEXE -p ./$SIMPR -f $PARFILE > screen_$SIMPR 2> PRtime
else
$TIMEEXE -p -o PRtime ./$SIMPR -f $PARFILE > screen_$SIMPR
fi
#
#
#
cd ..
