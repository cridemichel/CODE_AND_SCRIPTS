JOBID="138"
if [ "$2" == "" ]
then
NINI="260"
NFIN="0"
else
NINI="$1"
NFIN="$2"
fi
INICNF="startSUS.cnf"
TOTSTPS="50000"
BAKSTEPS="1000"
BAKSTEPSASCII="500000"
SAVESTPS="1000"
CT="cnftool "
MR="nohup mosrun -J$JOBID "
DELN="1"
if [ "$3" == "1" ]
then
if [ "$MR" != "" ] 
then
echo $JOBID > MOS_JOB_ID
else
rm MOS_JOB_ID
fi
fi
EHC="swell"
PF="ellipsoid_flex_mc.par"
SP="./set_params.sh"
WT="0.5"
N=$NINI
NPCNF=`cat $INICNF| awk -F : '{if ($1=="parnum") print $2}'` 
T="0.9128"
ACT="0.00226396"
BL="T-$T-z-$ACT"
STOP="0"
while [ $STOP -eq 0 ]
do 
echo "N=" $N " NFIN=" $NFIN
DELP=$[${NPCNF}-${N}]
echo "delp= " $DELP
if [ $NFIN -lt $NINI ]
then
NMIN=$[$[N]-$[DELN]]
NMAX=$N
else
NMAX=$[$[N]+$[DELN]]
NMIN=$N
fi
DN="N_${NMIN}_${NMAX}"
if [ ! -e $DN ]
then
mkdir $DN
fi
echo "NMIN=" $NMIN " NMAX=" $NMAX
cp $INICNF $DN/
# REORDER CONF AND REMOVE PARTICLES HERE
cp $PF $DN
cd $DN
cp ../set_params.py ../set_params.sh ../set_one_param.sh .
cp ../$CT .
rm -f COORD_TMP*
./$CT -rp $DELP $INICNF > _aaa_
cp _aaa_ $INICNF
EN="SWHE_$DN"
rm _aaa_
cp ../$EHC .
ln -sf ./$EHC $EN
NPAV=$NMIN
if [ "$NPAV" == "0" ] 
then
NPAV="1"
fi
NEXC=`echo 5.0*$NPAV/100.0|bc -l| awk '{printf("%d",$1)}'`
if [ "$NEXC" == "0" ] 
then
NEXC="1"
fi
echo "NPAV=" $NPAV " NEXC=" $NEXC
$SP $PF susnmin $NMIN susnmax $NMAX targetAccept 0.5 adjstepsMC 20000 zetaMC $ACT npav $NPAV nexc $NEXC temperat $T bakSteps $BAKSTEPS stepnum $TOTSTPS VSteps $SAVESTPS  bakStepsAscii $BAKSTEPSASCII
$MR ./$EN -fa $PF > screen &
sleep $WT
cd ..
if [ $NFIN -lt $NINI ]
then
N=$[$[N]-$[DELN]]
if [ $N -le $NFIN ]
then
STOP="1"
fi
else
N=$[$[N]+$[DELN]]
if [ $N -ge $NFIN ]
then
STOP="1"
fi
fi
done
