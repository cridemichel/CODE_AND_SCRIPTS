# $1 =pressure
# $2 = temperature 
# $3 = steps
if [ "$1" = "" ]
then
echo "Syntax: sim1statepnt <pressure> <temperature> <steps>"
exit 0
fi
PRESS="$1"
if [ "$2" = "" ]
then
echo "You have to supply the temperature"
exit 0
fi
TEMP="$2"
if [ "$3" = "" ]
then 
STEPS="10000"
else
STEPS="$3"
fi
PARFILE="ellipsoid_flex_mc.par"
DIRSIM="P-$PRESS"
if [ ! -e $DIRSIM ]
then 
mkdir $DIRSIM
fi
cp $PARFILE $DIRSIM
cd $DIRSIM
rm -f COORD_TMP*
rm -f Store-*
PRESS="$1"
ELLEXE="../ellipsHC"
INITEMP="2.0"
SIMPR="HCNPT-X0-P${PRESS}-T${TEMP}"
MOSRUN="mosrun"
#per ora il salvataggio è lineare
#=========== >>> PARAMETRI <<< =============
STORERATE="50.0"
USENNL=1
INIFILE="start.cnf"
PARNUM="1000"
#============================================
#N.B. it's supposed that we use NNL here!!
cp ../$INIFILE .
#==================================================================
echo "Simulating T=" $TEMP " P=" $PRESS
../set_params.sh $PARFILE temperat $TEMP ensembleMC 1 useNNL 0 rcut -1 P $PRESS bakStepsAscii 5000 stepnum $STEPS inifile $INIFILE endfile ${SIMPR}.cor
#../set_params.sh $PARFILE inifile start.cnf
ln -sf $ELLEXE $SIMPR
$MOSRUN ./$SIMPR -fa ./$PARFILE > screen_$SIMPR &
sleep 1
cd ..
