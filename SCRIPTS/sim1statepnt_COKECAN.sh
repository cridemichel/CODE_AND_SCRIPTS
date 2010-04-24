# $1 = volume fraction
# $2 = cicli di produzione (1 ciclo = 1 tempo di equilibratura) (se < 0 non fa la produzione)
# $3 = elongazione
# $4 = se < 0 non fa l'equilibratura, se = 0 usa il MSD per terminare l'equilibratura, se > 0 fa 
#      $4 passi di equilibratura ( 0 = default ) 
# $5 = temperatura
# $6 = sigma
if [ "$1" = "" ]
then
echo "Syntax: sim1statepnt <Volume_Fraction> <production_cycles> [elongation] [equilibration_steps] [temperature] [sigma]"
exit 0
fi
if [ "$2" = "" ]
then
echo "You have to supply the number of production cycles"
exit 0
fi
if [ "$3" = "" ]
then 
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
#echo "You have to supply an elongation!"
#echo "Elongazione: $EL"
#exit 0
else
EL=$3
fi
if [ "$4" == "" ]
then
EQSTPS=0
else
EQSTPS=$4
fi
if [ "$5" == "" ]
then
TEMP="0.15"
else
TEMP="$5"
fi
PARFILE=ellips.par
cp $PARFILE Phi$1
cd Phi$1
rm -f COORD_TMP*
rm -f Store-*
ELLEXE="../../ellipsoid"
if [ "$5" == "" ]
then
SIGMA="0.9"
else
SIGMA="$5"
fi
SIMRA="ell${EL}RA$1T${INITEMP}SIG$SIGMA"
SIMGR="ell${EL}GR$1T${INITEMP}SIG$SIGMA"
SIMEQ="ell${EL}EQ$1T${TEMP}SIG$SIGMA"
SIMPR="ell${EL}PR$1T${TEMP}SIG$SIGMA"
#per ora il salvataggio è lineare
#=========== >>> PARAMETRI <<< =============
STORERATE="50.0"
USENNL=1
INITEMP="2.0"
INIFILE="startIniSQ.cnf"
PARNUM=512
DT="0.05"
NNPR=1
BASEPR=1
INTSUM="20.0"
#============================================
#N.B. it's supposed that we use NNL here!!
PROL=`echo $EL | awk '{if ($0 >= 1.0) printf("1"); else printf("0");}'`
if [ $PROL -eq 1 ]
then
A0=$EL
B0=1.0
C0=1.0
RNNL=0.3
#INIL=`echo "2.0*e(1.0/3.0*l($PARNUM))*$A0" | bc -l`
#echo "qui INIL=" $INIL
if [ $USENNL -eq 0 ]
then
RCUT=`echo "2.0*$A0*1.01" | bc -l` 
fi
else
A0="1.0"
B0="$EL" 
C0="$EL"
RNNL=0.3
#INIL=`echo "5.0*e(1.0/3.0*l($PARNUM))*$B0" | bc -l`
if [ $USENNL -eq 0 ]
then
RCUT=`echo "2.0*$B0*1.01" | bc -l` 
NNLPAR="0"
fi
fi
if [ $USENNL -eq 1 ]
then
#usa le NNL con sticky spots!
NNLPAR="3"
RCUT=`echo "2.0*e(0.5*l(($A0+$RNNL)*($A0+$RNNL)+($B0+$RNNL)*($B0+$RNNL)+($C0+$RNNL)*($C0+$RNNL)))" | bc -l`
fi
echo "RCUT=" $RCUT " " "A=" $A0 "B=" $B0 "C=" $C0 "RNNL=" $RNNL "EL=" $EL
#RANDOMIZZAZIONE INIZIALE
cp $PARFILE rand_$PARFILE
#echo "L:" $INIL >> rand_$PARFILE
#>>> SET TEMPERATURE TO 1.0
../set_params.sh rand_$PARFILE stepnum 50 useNNL $NNLPAR temperat $INITEMP targetPhi 0.0 storerate 0.0 intervalSum 0.2 scalevel 1 rcut $RCUT rNebrShell $RNNL endfile ${SIMRA}.cor parnum $PARNUM inifile $INIFILE
ln -sf $ELLEXE $SIMRA
$SIMRA -f ./rand_${PARFILE} > screen_$SIMRA 
#exit 
#>>> EQUILIBRATE STARTING DENSITY (I.E. RANDOMIZE)
../set_params.sh rand_$PARFILE stepnum 200 targetPhi 0.0 storerate 0.0 intervalSum 2.0 scalevel 0 inifile ${SIMRA}.cor endfile ${SIMRA}.cor
ln -sf $ELLEXE $SIMRA 
$SIMRA -f ./rand_${PARFILE} >> screen_$SIMRA 
#CRESCITA
../set_params.sh $PARFILE stepnum 1000000 useNNL $NNLPAR targetPhi $1 storerate 0.0 intervalSum 0.05 rcut $RCUT rNebrShell $RNNL inifile ${SIMRA}.cor endfile ${SIMGR}.cor parnum $PARNUM 
ln -sf $ELLEXE $SIMGR
$SIMGR -f ./$PARFILE > screen_$SIMGR 
#exit 
#EQUILIBRATURA
if [ $EQSTPS -ge 0 ]
then
if [ $EQSTPS -eq 0 ]
then
STPS=10000000
TMSD=`echo "2.0*e((1.0/3.0)*l($A0*$B0*$C0))" | bc -l`
RMSD=-1.0 #1.57
else
STPS=$EQSTPS
TMSD="-1.0"
RMSD="-1.0"
fi
../set_params.sh $PARFILE stepnum $STPS targetPhi 0.0  temperat $TEMP storerate 0.0 intervalSum $INTSUM rmsd2end $RMSD tmsd2end $TMSD inifile ${SIMGR}.cor endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMEQ 
$SIMEQ -f ./$PARFILE > screen_$SIMEQ 
fi
#PRODUZIONE
if [ $2 -gt 0 ]
then
if [ $EQSTPS -ge 0 ]
then
if [ $EQSTPS -eq 0 ]
then
STCI=`cat screen_$SIMEQ | awk '{if ($1=="[MSDcheck]") print $3}'`
STPS=`echo "$STCI*$2"| bc -l`
else
STCI=$EQSTPS
STPS=`echo "$STCI*$2"| bc -l`
fi
else
STCI=$2
STPS=$2
fi
#NN=`echo "1+l($DT*$STCI/$STORERATE)/l(1.3)" | bc -l | awk '{printf("%d",$0)}'`
../set_params.sh $PARFILE stepnum $STPS targetPhi 0.0 base $BASEPR storerate $STORERATE intevalSum $INTSUM rmsd2end -1.0 tmsd2end -1.0 NN $NNPR inifile ${SIMEQ}.cor endfile ${SIMPR}.cor
ln -sf $ELLEXE $SIMPR
$SIMPR -f ./$PARFILE > screen_$SIMPR 
fi
cd ..
