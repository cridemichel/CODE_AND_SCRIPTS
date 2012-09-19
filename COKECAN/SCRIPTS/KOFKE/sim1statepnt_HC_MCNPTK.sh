# $1 =pressure
# $2 = temperature 
# $3 = steps
VSTEPS="500"
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
../set_params.sh $PARFILE temperat $TEMP ensembleMC 1 useNNL 0 rcut -1 P $PRESS bakStepsAscii 5000 stepnum $STEPS inifile $INIFILE endfile ${EXENAME}.cor VSteps $VSTEPS
#../set_params.sh $PARFILE inifile start.cnf
ln -sf $ELLEXE $EXENAME
$MOSRUN ./$EXENAME -fa ./$PARFILE > screen_$SIMPR &
if [ "$MOSRUN" = "mosrun" ] 
then
sleep 5
else
sleep 2
fi
