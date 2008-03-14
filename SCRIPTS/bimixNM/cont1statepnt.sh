PARFILE=BMsoft.par
if [ "$1" == "" ]
then
echo "con1statepnt.sh <temperature> <factor> <custom string>"
exit
fi
#cp $PARFILE T$1
cd T$1
rm -f COORD_TMP*
ELLEXE="../../bimix"
SIMCONT="bimixNM${3}T$1"
STORERATE="0.01"
PARNUM=1000
DT="0.002"
if [ "$2" == "" ]
then
FACT=5
else
FACT=$2
fi
LASTFILE=`ls Cnf*|sort -t - -k 3 -n | tail -1`
TOTSTPS=`cat $LASTFILE | awk -v fact=$FACT -F : '{if ($1=="totStep") print $2*fact}'`
echo "totSteps=>" $TOTSTPS
cat $LASTFILE | awk -v tstps=$TOTSTPS -F : '{if ($1=="totStep") print ("totStep:",tstps); else print $0}' > restart.cnf
#../set_params.sh $PARFILE stepnum $NSTPS inifile ${SIMGR}.cor endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMCONT 
$SIMCONT -ca  restart.cnf > screen_$SIMCONT
cd ..
