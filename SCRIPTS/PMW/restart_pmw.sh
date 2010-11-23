# $1 = temperature $2 = tau_alpha (in reduced units) $3 = simulation label $4 = dtTra $5=initial file
SETPARAMS="../../set_params.sh"
TEMP="$1"
PF="PMW.par"
LNPR="PMW-T$1-PR"
LNEQ="PMW-T$1-EQ"
EXEF="../../PMW"
INIF="restart.res"
TAUALPHA="$2"
PREEXE="/bin/mosrun"
#PREEXE=""
# equilibration time will be TAUALPHA*EQTA
EQTA="1"
# production run time will be TAUALPHA*PRTA
PRTA="1000"
STORERATEPR="0.01"
BASEPR="1.3"
NNPR=`echo $TAUALPHA $STORERATEPR $BASEPR | gawk '{printf("%d", 1+log($1/$2)/log($3));}'`
echo "BASEPR="  $BASEPR "NNPR=" $NNPR
#translational and rotational time step to fulfill DSE and SE
if [ $4 == "" ]
then
TRATS="0.212207"
else
AA=$4
if [ $[AA] -lt 0 ]
then
TRATS=`echo "0.212207 $AA"| awk '{print $1/(-$2)}'`
else
TRATS="$4"	
fi
fi
ROTTS=`echo "0.063662 0.212207 $TRATS" | awk '{print $3*$1/$2}'` 
#10 * tau_alpha
STEPS=`echo "" | gawk -v ta=$TAUALPHA -v trats=$TRATS -v eqta=$EQTA '{printf("%d",eqta*ta/trats);}'` 
echo "TEMP=$1 STEPS=$STEPS STORERATE=$STORERATEPR TAUALFA=$TAUALPHA" 
EQFN=T${TEMP}
if [ "$3" != "" ]
then	
FOLN=T${TEMP}-$3
else
FOLN=T${TEMP}-R	
fi
if [ ! -e $FOLN ]
then
mkdir $FOLN
fi
#choose Store file
if [ "$5" == "" ]
then
cd $EQFN
RF=`ls Store-*-*  | sort -t - -k 2 -k 3 -n | tail -1`
cp $RF ../$FOLN
cp $PF ../$FOLN
cd ..
else
INIF="$5"
cp $PF $FOLN
cp $INIF $FOLN
fi
cd $FOLN
#remove first header
if [ "$INIF" == "" ]
then
gunzip -f $RF
RF=`ls Store*`
cat $RF | awk 'BEGIN {cc=0} {if (cc==1) print $0; if ($1=="@@@") cc=1;}' > $INIF
rm $RF
fi
VS=`echo $TAUALPHA 100 $TRATS | awk '{printf("%d",$1/$2/$3)}'`
#PRODUZIONE
STEPS=`echo "" | gawk -v ta=$TAUALPHA -v trats=$TRATS -v prta=$PRTA '{printf("%d",prta*ta/trats);}'` 
$SETPARAMS $PF  intervalSum 50.0 storerate $STORERATEPR base $BASEPR NN $NNPR Dt $TRATS DtR $ROTTS stepnum  $STEPS inifile $INIF baksteps 0 VSteps $VS temperat $TEMP
ln -sf $EXEF $LNPR 
$PREEXE ./$LNPR -fa $PF > screenPR
cd ..
