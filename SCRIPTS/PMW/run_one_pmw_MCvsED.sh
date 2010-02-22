# $1 = temperature $2 = tau_alpha (in reduced units)
SETPARAMS="../set_params.sh"
TEMP="$1"
PF="PMW.par"
LN="PMW-T$1"
LNEQ="PMT-T$1-EQ"
EXEF="../../PMW"
INIF="restart.res"
TAUALPHA="$2"
# equilibration time will be TAUALPHA*EQTA
EQTA="2"
# production run time will be TAUALPHA*PRTA
PRTA="10"
STORERATEPR="0.01"
BASEPR="1.3"
NNPR=`echo $TAUALPHA $STORERATEPR $BASEPR | gawk '{printf("%d", 1+log($1/$2)/log($3);}'`
#translational and rotational time step to fulfill DSE and SE
TRATS="0.212207"
ROTTS="0.063662"
#10 * tau_alpha
STEPS=`echo "" | gawk -v ta=$TAUALPHA -v trats=$TRATS -v eqta=$EQTA '{printf("%d",eqta*ta/trats);}'` 
echo "TEMP=$1 STEPS=$STEPS STORERATE=$STORERATE TAUALFA=$TAUALPHA" 
if [ ! -e T$TEMP ]
then
mkdir T$TEMP
fi
cd T$TEMP
cp ../$PF .
VS=`echo $TAUALPHA 100 $TS | awk '{printf("%d",$1/$2/$3)}'`
#EQUILIBRATURA
$SETPARAMS $PF  intervalSum 50.0 storerate 0.0 Dt $TRATS DtR $ROTTS stepnum  $STEPS inifile $INIF NN 1 base 1 baksteps 0 VSteps $VS temperat $TEMP
ln -sf $EXEF $LNEQ
./$LN -fa $PF > screenEQ &
#store potential energy to later check equilibration
cp energy.dat energyEQ.dat
#store parameter file for equilibration
cp $PF EQ-$PF
#PRODUZIONE
STEPS=`echo "" | gawk -v ta=$TAUALPHA -v trats=$TRATS -v prta=$PRTA '{printf("%d",prta*ta/trats);}'` 
$SETPARAMS $PF  intervalSum 50.0 storerate $STORERATEPR base $BASEPR NN $NNPR Dt $TRATS DtR $ROTTS stepnum  $STEPS inifile $INIF baksteps 0 VSteps $VS temperat $TEMP
ln -sf $EXEF $LN 
./$LN -fa $PF > screen &
cd ..
