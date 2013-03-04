# $1 =pressure
# $2 = inverse temperature
# $3 = executable name 
# $4 = steps
# $5 = 1 (default) -> adjust volume and transl/rot displac to have 50% acceptance, 0 -> do not do that
VSTEPS="50"
TARGVOL="0.50"
TARGRT="0.50"
ADJSTEPS="10000"
RESETRT="500"
RESETVOL="500"
alias awk='LANG=C awk'
if [ "$1" = "" ]
then
echo "Syntax: sim1statepnt <pressure> <temperature> <steps>"
exit 0
fi
PRESS="$1"
if [ "$2" = "" ]
then
echo "You have to supply the inverse temperature (beta)"
exit 0
fi
BETA="$2"
if [ "$3" = "" ]
then
echo "You have to supply the executable name"
exit 0
fi
EXENAME="$3"
if [ "$4" = "" ]
then 
STEPS="10000"
else
STEPS="$4"
fi
if [ "$4 = "" ]
then
DOADJ="1"
else
DOADJ="$4"
fi
PARFILE="ellipsoid_flex_mc.par"
DIRSIM="P-$PRESS"
INIFILE="start.cnf"
cp ../$PARFILE .
rm -f COORD_TMP*
rm -f Store-*
rm -f CorFinal
PRESS="$1"
ELLEXE="../ellipsHC"
TEMP=`echo "1.0/$BETA"| bc -l`
SIMPR=`echo $PRESS $TEMP | awk '{printf("HCNPT-X0-P%.5f-T%.4f",$1,$2)}'`
MOSRUN=""
#per ora il salvataggio è lineare
#=========== >>> PARAMETRI <<< =============
STORERATE="50.0"
USENNL=1
PARNUM="1000"
#============================================
#N.B. it's supposed that we use NNL here!!
#cp ../$INIFILE .
#==================================================================
echo "[ " EXE:$EXENAME " ] Simulating BETA=" $BETA " P=" $PRESS
if [ "$DOADJ" == "1" ]
then
../set_params.sh $PARFILE temperat $TEMP adjstepsMC $ADJSTEPS targetAccept $TARGVOL targetAcceptVol $TARGRT ensembleMC 1 useNNL 0 rcut -1 P $PRESS bakStepsAscii 100000 stepnum $STEPS resetacceptVol $RESETVOL resetaccept $RESETRT inifile $INIFILE endfile ${EXENAME}.cor VSteps $VSTEPS
#../set_params.sh $PARFILE inifile start.cnf
ln -sf $ELLEXE $EXENAME
$MOSRUN ./$EXENAME -fa ./$PARFILE > screen_$SIMPR 
else
../set_params.sh $PARFILE temperat $TEMP targetAccept -1 targetAcceptVol -1 ensembleMC 1 useNNL 0 rcut -1 P $PRESS bakStepsAscii 100000 stepnum $STEPS inifile $INIFILE endfile ${EXENAME}.cor VSteps $VSTEPS
#../set_params.sh $PARFILE inifile start.cnf
ln -sf $ELLEXE $EXENAME
$MOSRUN ./$EXENAME -fa ./$PARFILE > screen_$SIMPR &
fi
if [ "$MOSRUN" = "mosrun" ] 
then
sleep 5
else
sleep 5
fi
