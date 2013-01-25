JOBID="138"
NINI="1550"
NFIN="1400"
INICNF="startSUS.cnf"
CT="cnftool"
MR="mosrun -J$JOBID "
if [ "$MR" != "" ] 
then
echo $JOBID > MOS_JOB_ID
else
rm MOS_JOB_ID
fi
EHC="ellipsHC"
PF="ellipsoid_flex_mc.par"
SP="../set_params.sh"
WT="2"
N=$NINI
NPCNF=`cat $INICNF| awk -F : '{if ($1=="parnum") print $2}'` 
T="0.10"
ACT="0.009"
BL="T-$T-z-$ACT"
while [ $N -gt $NFIN ]
do 
echo "N=" $N " NFIN=" $NFIN
DELP=$[${NPCNF}-${N}]
echo "delp= " $DELP
NMIN=$[$[N]-1]
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
../$CT -rc 2  -rp $DELP $INICNF > _aaa_
cp _aaa_ $INICNF
EN="HC_$DN"
rm _aaa_
ln -sf ../$EHC $EN
NPAV=$NMIN
NEXC=`echo 5.0*$NPAV/100.0|bc -l| awk '{printf("%d",$1)}'`
echo "NPAV=" $NPAV " NEXC=" $NEXC
$SP $PF susnmin $NMIN susnmax $NMAX zetaMC $ACT npav $NPAV nexc $NEXC
$MR ./$EN -fa $PF > screen &
sleep $WT
cd ..
N=$[${N}-1]
done
