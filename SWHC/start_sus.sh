JOBID="138"
NINI="500"
NFIN="0"
INICNF="startSUS.cnf"
BAKSTEPS="5000"
CT="cnftool "
MR="nohup mosrun -J$JOBID "
DELN="1"
if [ "$MR" != "" ]
then
  echo $JOBID > MOS_JOB_ID
else
rm MOS_JOB_ID
fi
EHC="SWHC"
PF="ellipsoid_flex_mc.par"
SP="../set_params.sh"
WT="1"
N=$NINI
NPCNF=`cat $INICNF| awk -F : '{if ($1=="parnum") print $2}'` 
T="0.7"
ACT="1.05"
BL="T-$T-z-$ACT"
while [ $N -gt $NFIN ]
do 
echo "N=" $N " NFIN=" $NFIN
DELP=$[${NPCNF}-${N}]
echo "delp= " $DELP
NMIN=$[$[N]-$[DELN]]
NMAX=$N
DN="N_${NMIN}_${NMAX}"
if [ ! -e $DN ]
then
mkdir $DN
fi
cp $INICNF $DN/
# REORDER CONF AND REMOVE PARTICLES HERE
cp $PF $DN
cd $DN
rm -f COORD_TMP*
../$CT -rp $DELP $INICNF > _aaa_
cp _aaa_ $INICNF
EN="HC_$DN"
rm _aaa_
ln -sf ../$EHC $EN
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
$SP $PF susnmin $NMIN susnmax $NMAX zetaMC $ACT npav $NPAV nexc $NEXC temperat $T bakSteps $BAKSTEPS
$MR ./$EN -fa $PF > screen &
sleep $WT
cd ..
N=$[${N}-${DELN}]
done
