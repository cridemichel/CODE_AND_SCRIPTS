# $1 = temperature $2 = tau_alpha (in reduced units)
SETPARAMS="../../set_params.sh"
TEMP="$1"
PF="PMW.par"
LNPR="PMW-T$1-PR"
LNEQ="PMT-T$1-EQ"
EXEF="../../PMW"
INIF="restart.res"
TAUALPHA="$2"
# equilibration time will be TAUALPHA*EQTA
EQTA="1"
# production run time will be TAUALPHA*PRTA
PRTA="10"
STORERATEPR="0.01"
BASEPR="1.3"
NNPR=`echo $TAUALPHA $STORERATEPR $BASEPR | gawk '{printf("%d", 1+log($1/$2)/log($3));}'`
echo "BASEPR="  $BASEPR "NNPR=" $NNPR
#translational and rotational time step to fulfill DSE and SE
TRATS="0.212207"
ROTTS="0.063662"
#10 * tau_alpha
STEPS=`echo "" | gawk -v ta=$TAUALPHA -v trats=$TRATS -v eqta=$EQTA '{printf("%d",eqta*ta/trats);}'` 
echo "TEMP=$1 STEPS=$STEPS STORERATE=$STORERATEPR TAUALFA=$TAUALPHA" 
if [ ! -e T$TEMP ]
then
mkdir T$TEMP
fi
cd T$TEMP
cp ../$INIF .
cp ../$PF .
VS=`echo $TAUALPHA 100 $TRATS | awk '{printf("%d",$1/$2/$3)}'`
#EQUILIBRATURA
$SETPARAMS $PF  intervalSum 50.0 storerate -1.0 Dt $TRATS DtR $ROTTS stepnum  $STEPS inifile $INIF NN 1 base 1 baksteps 0 VSteps $VS temperat $TEMP
ln -sf $EXEF $LNEQ
./$LNEQ -fa $PF > screenEQ 
#store potential energy to later check equilibration
cp energy.dat energyEQ.dat
#store parameter file for equilibration
cp $PF EQ-$PF
#PRODUZIONE
STEPS=`echo "" | gawk -v ta=$TAUALPHA -v trats=$TRATS -v prta=$PRTA '{printf("%d",prta*ta/trats);}'` 
$SETPARAMS $PF  intervalSum 50.0 storerate $STORERATEPR base $BASEPR NN $NNPR Dt $TRATS DtR $ROTTS stepnum  $STEPS inifile $INIF baksteps 0 VSteps $VS temperat $TEMP
ln -sf $EXEF $LNPR 
./$LNPR -fa $PF > screenPR
cd ..
