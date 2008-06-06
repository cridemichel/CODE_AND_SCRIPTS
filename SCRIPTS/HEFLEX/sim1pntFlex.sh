if [ "$1" = "" ]
then
echo "Syntax: sim1pntFlex [equilibration_steps]"
exit 0
fi
if [ "$1" = "" ]
then
EQSTPS=4000
else
EQSTPS=$1
fi
N2=`ls -1 conf-1-*.dat|wc -l`
NT=`ls -1 conf-*-*.dat | wc -l`
N1=`echo $NT/$N2 | bc`
n1="0"
n2="0"
#densità superficiale di antigeni
CDIR=`pwd`
BN=`basename $CDIR`
SIG=`echo $BN | awk -F 'sigma' '{print $2}'`
N1MAX="5"
N2MAX="5"
if [ $N1 -gt $N1MAX ]
then
N1=$N1MAX
fi
if [ $N2 -gt $N2MAX ]
then
N2=$N2MAX
fi
echo "N1=" $N1 "N2=" $N2
while [ $n1 -lt $N1 ]
do
n2="0"
while [ $n2 -lt $N2 ]
do
#echo "n1=" $n1 "n2=$n2"
if [ ! -d RUN-${n1}-${n2} ]
then
mkdir RUN-${n1}-${n2}
fi
cd RUN-${n1}-${n2}
PARFILE="ellipsoid_flex.par"
INICONF="conf-$n1-$n2.dat"
rm -f COORD_TMP*
rm -f Store-*
ELLEXE="../../ellipsoid"
SIMEQ="ellflexEQ-${SIG}-${n1}-${n2}"
SIMPR="ellflexPR-${SIG}-${n1}-${n2}"
STORERATE="1000.0"
USENNL=0
PARNUM=2024
DT="0.05"
RCUT="-1"
#N.B. it's supposed that we use NNL here!!
cp ../$PARFILE .
cp ../$INICONF .
#elimina l'interazione antigene-anticorpo durante l'equilibratura
cat $INICONF | awk 'BEGIN { nat=0 } {if (!((($1=="0" && $2=="1")||($1=="1" && $2=="1")) && (nat == 2))) print $0; if ($0=="@@@") nat++;}' > iniconf.dat
if [ $EQSTPS -eq 0 ]
then
STPS=4000
else
STPS=$EQSTPS
fi
../set_params.sh $PARFILE stepnum $STPS storerate 0.0 scalevel 1 rescaleTime 5.0 intervalSum 10.0 inifile iniconf.dat endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMEQ 
$SIMEQ -fa ./$PARFILE > screen_$SIMEQ 
#PRODUZIONE
I1="0 1 5 0 1 0 100000 1"
I2="1 1 5 0 1 0 100000 1"
#attiva l'interazione antigene-anticorpo
cat CorFinal | awk -v i1="$I1" -v i2="$I2" 'BEGIN {nat=0} {if ($0=="@@@") nat++; if ($0=="@@@" && nat==3) {print i1; print i2; print "@@@"} else print $0}' > CorIni
../set_params.sh $PARFILE stepnum 10000000 base 1 NN 1 storerate $STORERATE intervalSum 10.0 inifile CorIni endfile ${SIMPR}.cor
ln -sf $ELLEXE $SIMPR
$SIMPR -fa ./$PARFILE > screen_$SIMPR 
cd ..
n2=$[$n2+1]
done
n1=$[$n1+1]
done
