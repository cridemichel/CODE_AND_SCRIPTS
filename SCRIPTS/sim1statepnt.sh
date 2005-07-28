# $1 = volume fraction
# $2 = cicli di produzione (1 ciclo = 1 tempo di equilibratura) (se < 0 non fa la produzione)
# $3 = elongazione
# $4 = se < 0 non fa l'equilibratura, se = 0 usa il MSD per terminare l'equilibratura, se > 0 fa 
#      $4 passi di equilibratura ( 0 = default ) 
if [ $3 = "" ]
then 
echo "You have to supply an elongation!"
exit 0
fi
if [ $4 = "" ]
then
EQSTPS=0
fi
EL=$3
PARFILE=ellips.par
cp $PARFILE Phi$1
cd Phi$1
ELLEXE="../../ellipsoid"
SIMRA="ell${EL}RA$1"
SIMGR="ell${EL}GR$1"
SIMEQ="ell${EL}EQ$1"
SIMPR="ell${EL}PR$1"
STORERATE="0.01"
USENNL=1
#N.B. it's supposed that we use NNL here!!
if [ $EL -gt 1.0 ]
then
A0=$EL
B0=1.0
C0=1.0
RNNL=0.12
if [ $USENNL -eq 0 ]
then
RCUT=`echo "2.0*$A0*1.01" | bc -l` 
fi
else
A0=1.0
B0=`echo "1.0/$EL" | bc -l`
C0=`echo "1.0/$EL" | bc -l`
RNNL=0.2
if [ $USENNL -eq 0 ]
then
RCUT=`echo "2.0*$B0*1.01" | bc -l` 
fi
fi
if [ $USENNL -eq 1 ]
then
RCUT=`echo "2.0*e(0.5*l(($A0+$RNNL)*($A0+$RNNL)+($B0+$RNNL)*($B0+$RNNL)+($C0+$RNNL)*($C0+$RNNL)))" | bc -l`
fi
#RANDOMIZZAZIONE INIZIALE
cp $PARFILE rand_$PARFILE
INIL=`echo "" | bc -l`
echo "L:" $INIL >> rand_$PARFILE
#>>> SET TEMPERATURE TO 1.0
../set_params.sh rand_$PARFILE stepnum 10 targetPhi 0.0 storerate 0.0 intervalSum 0.2 scalevel 1 rcut $RCUT rNebrShell $RNNL inifile '*' endfile ${SIMRA}.cor
ln -sf $ELLEXE $SIMRA
$SIMRA -f ./rand_${PARFILE} > screen_$SIMRA 
#>>> EQUILIBRATE STARTING DENSITY (I.E. RANDOMIZE)
../set_params.sh rand_$PARFILE stepnum 1000 targetPhi 0.0 storerate 0.0 intervalSum 0.2 scalevel 0 
ln -sf $ELLEXE $SIMRA inifile ${SIMRA}.cor endfile ${SIMRA}.cor
$SIMRA -f ./rand_${PARFILE} >> screen_$SIMRA 
#CRESCITA
../set_params.sh $PARFILE stepnum 1000000 targetPhi $1 storerate 0.0 intervalSum 0.02 rcut $RCUT rNebrShell $RNNL inifile ${SIMRA}.cor endfile ${SIMGR}.cor
ln -sf $ELLEXE $SIMGR
$SIMGR -f ./$PARFILE > screen_$SIMGR 
#EQUILIBRATURA
if [ $EQSTPS -gt 0 ]
then
if [ $EQSTPS -eq 0 ]
STPS=10000000
TMSD=`echo "2.0*e((1.0/3.0)*l($A0*$B0*$C0))" | bc -l`
RMSD=3.14
else
STPS=$EQSTPS
TMSD="-1.0"
RMSD="-1.0"
fi
../set_params.sh $PARFILE stepnum $STPS targetPhi 0.0  storerate 0.0 intervalSum 2.0 rmsd2end $RMSD tmsd2end $TMSD inifile ${SIMGR}.cor endfile ${SIMEQ}.cor
ln -sf $ELLEXE $SIMEQ 
$SIMEQ -f ./$PARFILE > screen_$SIMEQ 
fi
#PRODUZIONE
if [ $2 -gt 0 ]
then
STCI=`cat screen_$SIMEQ | awk '{if ($1=="[MSDcheck]") print $3}'`
STPS=`echo "$STCI*$2"| bc -l`
NN=`echo "l($STCI/$STORERATE)/l(10.0)" | bc -l | awk 'printf("%d",$0)'`
../set_params.sh $PARFILE stepnum $STPS targetPhi 0.0 storerate $STORERATE intevalSum 2.0 rmsd2end -1.0 tmsd2end -1.0 NN $NN inifile ${SIMEQ}.cor endfile ${SIMPR}.cor
ln -sf $ELLEXE $SIMPR
$SIMPR -f ./$PARFILE > screen_$SIMPR 
fi
cd ..
