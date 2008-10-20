PARFILE=BMsoft.par
if [ "$5" == "" ]
then
echo "You must supply the temperature, the number of cycles and optionally a custom string, i.e.:"
echo "sim1statepnt.sh <temperature> <num. of cycles> <custom string> <Sim. Num.> <eq. steps>"
exit
fi
ln -sf $HOME/MDsimul/bin/bimix bimixNM 
if [ ! -d T$1-$4 ]
then
mkdir T$1-$4
fi
cp $PARFILE T$1-$4/
cd T$1
rm -f COORD_TMP*
BMEXE="../bimixNM"
SIMRA="bimixNM$3RA$1"
SIMEQ="bimixNM$3EQ$1"
SIMPR="bimixNM$3PR$1"
STORERATE="0.01"
SIG2EQ="2.0"
USENNL=1
PARNUM=1000
DT="0.002"
#RANDOMIZZAZIONE INIZIALE
cp $PARFILE rand_$PARFILE
#>>> SET TEMPERATURE TO 10.0
../set_params.sh rand_$PARFILE endFormat 2 Nose 2 steplength $DT stepnum 1000 chkeqstps 50 eqFact 5.0 endfile ${SIMRA}.cor parnum $PARNUM temperat 10.0 sResetSteps 100 CMreset 100 bakSaveMode 0 bakStepsAscii 500000
ln -sf $BMEXE $SIMRA
$SIMRA -f ./rand_${PARFILE} > screen_$SIMRA 
cp rand_$PARFILE $PARFILE
#EQUILIBRATURA
../set_params.sh $PARFILE Nose 2 stepnum $5 inifile ${SIMRA}.cor endfile ${SIMEQ}.cor sResetSteps 1000 CMreset 0 chkeqstps 50 eqFact $SIG2EQ temperat $1
ln -sf $BMEXE $SIMEQ 
$SIMEQ -f ./$PARFILE > screen_$SIMEQ 
#PRODUZIONE
if [ ! -d ../T$1-$4 ]
then
mkdir ../NVE-T$1-$4
fi
cp screen_$SIMEQ ../NVE-T$1-$4
cp $SIMEQ.cor ../NVE-T$1-$4
cp BMsoft.par ../NVE-T$1-$4
cd ../NVE-T$1-$4
rm -f COORD_TMP*
rm -f Cnf*
#STCI=`cat screen_$SIMEQ | awk '{if ($1=="[MSDcheck]") print $3}'`
STPS=`echo "$5*$2"| bc -l`
NN=`echo "l($STCI)/l(1.4)" | bc -l | awk '{printf("%d",$0)}'`
../set_params.sh $PARFILE stepnum $STPS chkeqstps 0 NN $NN inifile ${SIMEQ}.cor endfile ${SIMPR}.cor sResetSteps 0 bakSaveMode 0
ln -sf $BMEXE $SIMPR
./$SIMPR -f ./$PARFILE > screen_$SIMPR  
cd ..
