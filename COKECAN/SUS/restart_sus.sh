JOBID="138"
NINI="2050"
NFIN="1750"
INICNF="startSUS.cnf"
CT="cnftool"
DELN="1"
BAKSTEPS="5000"
MR="nohup mosrun -J$JOBID "
if [ "$MR" != "" ] 
then
echo $JOBID > MOS_JOB_ID
else
rm MOS_JOB_ID
fi
EHC="ellipsHC"
PF="ellipsoid_flex_mc.par"
SP="../set_params.sh"
WT="1"
N=$NINI
NPCNF=`cat $INICNF| awk -F : '{if ($1=="parnum") print $2}'` 
T="0.135"
ACT="1.0708"
BL="T-$T-z-$ACT"
while [ $N -gt $NFIN ]
do 
echo "N=" $N " NFIN=" $NFIN
NMIN=$[$[N]-$[DELN]]
NMAX=$N
DN="N_${NMIN}_${NMAX}"
cd $DN
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
N=$[${N}-${DELN}]
continue
fi
rm -f Cnf*
rm -f COORD_TMP[0,1]
EN="HC_$DN"
#ln -sf ../$EHC $EN
#NPAV=$NMIN
#NEXC=`echo 5.0*$NPAV/100.0|bc -l| awk '{printf("%d",$1)}'`
#echo "NPAV=" $NPAV " NEXC=" $NEXC
$SP $PF zetaMC $ACT flip_prob -0.1 bakSteps $BAKSTEPS
#ADD FLIP MOVE
#echo "flip_prob: 0.1" >> $PF
$MR ./$EN -fa $PF > screen &
sleep $WT
cd ..
N=$[${N}-${DELN}]
done
