# $1 = initial configuration
# $2 = total steps 
# $3 = equilibration steps
# $4 = rcut
# $5 = extra label
# $6 = volume fraction
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
RCUT="-1"
else
RCUT="$4"
fi
if [ "$5" = "" ]
then
EXTLAB="SE"
else
EXTLAB="$5"
fi
INIFILE="$1"
TEMP="1.0"
TIMEEXE="/usr/bin/time"
PARFILE=ellipsoid_flex.par
#get elongation from directory name
EL=`pwd | awk -F 'X0_' '{print $2}'`
if [ "$6" == "" ]
then
PHI=`echo $INIFILE | awk -F 'Phi' '{print $2}'| awk -F _ '{print $1}'`
else
PHI="$6"
fi
SETPARAMS="../../set_params.sh"
PD="PHI_${PHI}"
if [ ! -e $PD ]
then
mkdir $PD
fi
cp $PARFILE $PD
cp $INIFILE $PD
#exit
cd $PD
rm -f COORD_TMP*
rm -f Store-*
ELLEXE="../ellipsoid"
SIMPR="supell_X0_${EL}_Phi${PHI}_${EXTLAB}_PR"
SIMEQ="supell_X0_${EL}_Phi${PHI}_${EXTLAB}_EQ"
STORERATE="0.0"
USENNL=1
PARNUM=512
PARNUMA=512
DT="0.05"
RNNL="0.12"
#
#
#equilibration run
$SETPARAMS $PARFILE Dt $DT stepnum $EQSTP VSteps 0 temperat $TEMP scalevel 1 rescaleTime 0.5 targetPhi 0.0 storerate $STORERATE intervalSum 5.0 DtrCalc 0 rcut $RCUT rotMSDCalc 0 rmsd2end -1.0 tmsd2end -1.0 inifile $INIFILE endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMPR
./$SIMEQ -fa $PARFILE > screen_$SIMEQ 
#
#
#production run
$SETPARAMS $PARFILE stepnum $TOTSTP VSteps 0 temperat $TEMP targetPhi 0.0 storerate $STORERATE intervalSum 5.0 DtrCalc 0 rcut $RCUT rotMSDCalc 0 rmsd2end -1.0 tmsd2end -1.0 inifile ${SIMEQ}.cor endfile ${SIMPR}.cor
ln -sf $ELLEXE $SIMPR
$TIMEEXE -p -o PRtime ./$SIMPR -fa $PARFILE > screen_$SIMPR 
#
#
#
cd ..
