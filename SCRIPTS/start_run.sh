# $1 = volume fraction
# $2 = passi per equilibratura
# $3 = passi per produzione
# $4 = file di parametri
EL=3.0
PARFILE=ellips.par
cp $PARFILE Phi$1
cd Phi$1
ELLEXE="../ellipsoid230605"
SIMGR="ell${EL}GR$1"
SIMEQ="ell${EL}EQ$1"
SIMPR="ell${EL}PR$1"
#CRESCITA
../set_params.sh $PARFILE stepnum 50000 targetPhi $1 storerate 0.0 intervalSum 0.02 
ln -sf $ELLEXE $SIMGR
$SIMGR -f ./$PARFILE > screen_$SIMGR 
#EQUILIBRATURA
../set_params.sh $PARFILE stepnum $2 targetPhi 0.0  storerate 0.0
ln -sf $ELLEXE $SIMEQ 
$SIMEQ -f ./$PARFILE > screen_$SIMEQ 
#PRODUZIONE
../set_params.sh $PARFILE stepnum $3 targetPhi 0.0 storerate 0.01
ln -sf $ELLEXE $SIMPR
$SIMPR -f ./$PARFILE > screen_$SIMPR 
cd ..
