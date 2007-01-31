if [ "$1" = "" ]
then
echo "Syntax: sim1pntFlex [equilibration_steps]"
exit 0
fi
if [ "$1" = "" ]
then
EQSTPS=2000
else
EQSTPS=$1
fi
N2=`ls -1 CONF-1-*.dat|wc -l`
NT=`ls -1 CONF-*-*.dat | wc -l`
N1=`echo NT/N2 | bc`
nr="0"
#densità superificiale di antigeni
SIG=`basename `pwd`| awk -F 'sigma' '{print $2}'`
while [ $n1 -lt $N1 ]
do
while [ $n2 -lt $N2 ]
do
if [ ! -d RUN-${n1}-${n2} ]
then
mkdir RUN-${n1}-${n2}
fi
cd RUN-${n1}-${n2}
PARFILE="ellips_flex.par"
INICONF="CONF-$n1-$n2.dat"
rm -f COORD_TMP*
rm -f Store-*
ELLEXE="../../ellipsoid"
SIMEQ="ellflexEQ-${SIG}-${n1}-${n2}"
SIMPR="ellflexPR-${SIG}-${n1}-${n2}"
STORERATE="0.01"
USENNL=0
PARNUM=2024
DT="0.05"
RCUT="-1"
#N.B. it's supposed that we use NNL here!!
cp ../$PARFILE .
mv ../$INICONF .
if [ $EQSTPS -eq 0 ]
then
STPS=2000
else
STPS=$EQSTPS
fi
../../set_params.sh $PARFILE stepnum $STPS storerate 0.0 scalevel 1 rescaleTime 5.0 intervalSum 10.0 inifile $INICONF endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMEQ 
$SIMEQ -fa ./$PARFILE > screen_$SIMEQ 
#PRODUZIONE
../../set_params.sh $PARFILE stepnum $STPS storerate $STORERATE intervalSum 10.0 inifile ${SIMEQ}.cor endfile ${SIMPR}.cor
ln -sf $ELLEXE $SIMPR
$SIMPR -f ./$PARFILE > screen_$SIMPR 
cd ..
n2=$[$n2+1]
done
n1=$[$n1+1]
done
