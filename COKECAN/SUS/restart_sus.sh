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
CT="cnftool"
DELN="1"
BAKSTEPS="20000"
TOTSTPS="20000"
SAVESTPS="2000"
BAKSTEPSASCII="500000"
MR="nohup mosrun -J$JOBID "
if [ "$3" == "1" ] 
then
if [ "$MR" != "" ] 
then
echo $JOBID > MOS_JOB_ID
else
rm MOS_JOB_ID
fi
fi
EHC="ellipsHC"
PF="ellipsoid_flex_mc.par"
SP="./set_params.sh"
WT="0.1"
N=$NINI
NPCNF=`cat $INICNF| awk -F : '{if ($1=="parnum") print $2}'` 
T="0.548"
ACT="0.000919"
BL="T-$T-z-$ACT"
STOP="0"
while [ $STOP -eq 0 ]
do 
echo "N=" $N " NFIN=" $NFIN
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
echo "DIR= " $DN " does not exist, skipping..."
if [ $NFIN -lt $NINI ]
then
N=$[${N}-${DELN}]
else
N=$[${N}+${DELN}]
fi
continue
fi
cd $DN
cp ../set_params.py ../set_params.sh ../set_one_param.sh .
# FIND THE LATEST COORD_TMP_ASCII FILE (and check that it is not corrupted)
CFO=`ls -rt COORD_TMP_ASCII*| head -1`
CFN=`ls -rt COORD_TMP_ASCII*| tail -1`
PNO=`cat $CFO| awk -F : '{if ($1=="parnum") {print $2; exit;}}'`
PNN=`cat $CFN| awk -F : '{if ($1=="parnum") {print $2; exit;}}'`
POSATO=`cat $CFO | awk 'BEGIN{at=0} {if ($0=="@@@") {at++;}; if (at==3) {print NR; exit;}}'`
POSATN=`cat $CFN | awk 'BEGIN{at=0} {if ($0=="@@@") {at++;}; if (at==3) {print NR; exit;}}'`
WCO=`cat $CFO| wc -l`
WCN=`cat $CFN| wc -l`
EXPLO=`echo $POSATO+$PNO+1|bc`
EXPLN=`echo $POSATN+$PNN+1|bc`
echo "WCO=" $WCO " CFO=" $CFO " PNO=" $PNO " POSATO=" $POSATO " EXPLO=" $EXPLO
echo "WCN=" $WCN " CFN=" $CFN " PNN=" $PNN " POSATN=" $POSATN " EXPLN=" $EXPLN
echo "NMIN=" $NMIN " NMAX=" $NMAX
if [ $[EXPLN] -eq $[WCN] ]
then
echo "I chose CFN=" $CFN
cat $CFN | awk 'BEGIN {att=0} {if (att==1) print $0; if ($0=="@@@") att=1;}' > $INICNF
elif [ $[EXPLO] -eq $[WCO] ]
then
echo "I chose CFO=" $CFO
cat $CFO | awk 'BEGIN {att=0} {if (att==1) print $0; if ($0=="@@@") att=1;}' > $INICNF
else
echo "DIR=" $DN
echo "FATAL ERROR: both COORD_TMP_ASCII corrupted!"
echo "I skip this directory..."
cd ..
if [ $NFIN -lt $NINI ]
then
N=$[${N}-${DELN}]
else
N=$[${N}+${DELN}]
fi
continue
fi
rm -f Cnf*
rm -f COORD_TMP[0,1]
EN="HC_$DN"
#ln -sf ../$EHC $EN
#NPAV=$NMIN
#NEXC=`echo 5.0*$NPAV/100.0|bc -l| awk '{printf("%d",$1)}'`
#echo "NPAV=" $NPAV " NEXC=" $NEXC
$SP $PF zetaMC $ACT flip_prob -0.1 bakSteps $BAKSTEPS VSteps $SAVESTPS temperat $T stepnum $TOTSTPS bakStepsAscii $BAKSTEPSASCII
#ADD FLIP MOVE
#echo "flip_prob: 0.1" >> $PF
$MR ./$EN -fa $PF > screen &
sleep $WT
cd ..
if [ $NFIN -lt $NINI ]
then
N=$[${N}-${DELN}]
if [ $N -le $NFIN ]
then
STOP="1"
fi
else
N=$[${N}+${DELN}]
if [ $N -ge $NFIN ]
then
STOP="1"
fi
fi
done
