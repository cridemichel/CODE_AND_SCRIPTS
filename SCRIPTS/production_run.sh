# $1 = volume fraction
# $2 = cicli di produzione (1 ciclo = 1 tempo di equilibratura) (se < 0 non fa la produzione)
# $3 = elongazione
# $4 = se < 0 non fa l'equilibratura, se = 0 usa il MSD per terminare l'equilibratura, se > 0 fa 
#      $4 passi di equilibratura ( 0 = default ) 
if [ "$1" = "" ]
then
echo "Syntax: production_run.sh <Volume_Fraction> <production_cycles> <steps for each cycle> [elongation]"
exit 0
fi
if [ "$2" = "" ]
then
echo "You have to supply the number of production cycles"
exit 0
fi
if [ "$4" = "" ]
then 
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
#echo "You have to supply an elongation!"
#echo "Elongazione: $EL"
#exit 0
else
EL=$4
fi
if [ "$3" = "" ]
then
echo "you have to supply the steps of each cycle"
else
EQSTPS=$3
fi
PARFILE=ellipsPR.par
cp $PARFILE Phi$1-PR
cd Phi$1-PR
rm -f COORD_TMP*
ELLEXE="../../ellipsoid"
SIMPR="ell${EL}PR$1"
STORERATE="0.1"
NN=`echo `
USENNL=1
PARNUM=512
PARNUMA=512
A0=`cat ellips.res | awk -F ':' '{if ($1=="a") print $2}' | awk '{print $1}'`
B0=`cat ellips.res | awk -F ':' '{if ($1=="b") print $2}' | awk '{print $1}'`
C0=`cat ellips.res | awk -F ':' '{if ($1=="c") print $2}' | awk '{print $1}'`
DT="0.05"
#N.B. it's supposed that we use NNL here!!
PROL=`echo $EL | awk '{if ($0 >= 1.0) printf("1"); else printf("0");}'`
if [ $PROL -eq 1 ]
then
RNNL=0.15
#echo "qui INIL=" $INIL
if [ $USENNL -eq 0 ]
then
RCUT=`echo "2.0*$A0*1.01" | bc -l` 
fi
else
RNNL=0.15
if [ $USENNL -eq 0 ]
then
RCUT=`echo "2.0*$B0*1.01" | bc -l` 
fi
fi
if [ $USENNL -eq 1 ]
then
RCUT=`echo "2.0*e(0.5*l(($A0+$RNNL)*($A0+$RNNL)+($B0+$RNNL)*($B0+$RNNL)+($C0+$RNNL)*($C0+$RNNL)))" | bc -l`
fi
echo "RCUT=" $RCUT " A0=" $A0 " B0=" $B0 " C0=" $C0 " RNNL=" $RNNL "EL=" $EL
#RANDOMIZZAZIONE INIZIALE
#>>> SET TEMPERATURE TO 1.0
STCI=$EQSTPS
STPS=$EQSTPS
NN=`echo "1+l($DT*$STCI/$STORERATE)/l(1.3)" | bc -l | awk '{printf("%d",$0)}'`
../set_params.sh $PARFILE stepnum $STPS targetPhi 0.0 storerate $STORERATE intevalSum 5.0 rmsd2end -1.0 tmsd2end -1.0 NN $NN inifile ${SIMEQ}.cor endfile ${SIMPR}.cor inifile ellips.res rcut $RCUT
ln -sf $ELLEXE $SIMPR
$SIMPR -fa ./$PARFILE > screen_$SIMPR 
cd ..
