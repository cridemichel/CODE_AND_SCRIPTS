PARFILE=BMsoft.par
if [ "$1" == "" ]
then
echo "Syntax: NVEfromBVT.sh <temperature> [FACTOR] [NN] [custom string] " 
exit
fi
if [ ! -e NVE-T$1 ]
then
mkdir NVE-T$1
fi
cp $PARFILE NVE-T$1/
cd NVE-T$1
rm -f COORD_TMP*
rm -f Cnf*
ELLEXE="../../bimix"
SIMNVE="bimixNM-NVE-${4}T$1"
STORERATE="0.01"
PARNUM=1000
DT="0.002"
if [ "$2" == "" ]
then
FACT=1
else
FACT=$2
fi
LASTFILE=`ls ../T$1/Cnf*|sort -t - -k 3 -n | tail -1`
if `zcat $LASTFILE > /dev/null 2>&1` 
then
CATEXE="zcat"
else
CATEXE="cat"
fi
TOTSTPS=`$CATEXE $LASTFILE | awk -v fact=$FACT -F : '{if ($1=="totStep") print $2*fact}'`
$CATEXE $LASTFILE | awk -v tstps=$TOTSTPS -F : ' BEGIN {doout=0} {if (doout==1) print $0; if ($0=="@@@") doout=1;}' > restart.cor
echo "totSteps=>" $TOTSTPS
../set_params.sh $PARFILE bakSaveMode 1 stepnum $TOTSTPS inifile restart.cor endFormat 2 endfile ${SIMNVE}.cor Nose 0 nRun NVE
if [ "$3" != "" ]
then 
#NOTA: se si specifica NN FACT si riferisce ai blocchi di NN configurazioni 
#  invece se NN non viene specificato allora FACT è un fattore rispetto ai passi 
#  che gia' si trovano in BMsoft.par
TOTSTPS=`echo "(1.4^$3)*$FACT" | bc -l | awk '{printf("%d",$0)}'`
../set_params.sh $PARFILE NN $3 stepnum $TOTSTPS
echo "===> NN= " $3 " Total Steps=" $TOTSTPS
fi
ln -sf $ELLEXE $SIMNVE 
$SIMNVE -fa $PARFILE > screen_$SIMNVE
cd ..
