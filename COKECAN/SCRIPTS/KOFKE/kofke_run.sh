#To start Kofke run create the following files:
# BetaEnd.dat : 1/(final temparature)
# si deve anche creare il file $KFN con dentro i valori per gli stati iniziali, ossia:
# <beta> <pressure> <av. volume ISO> <av. energy ISO> <av. vol. NEM> <av. energy NEM> 
#
# $1 è il passo d'integrazione in 1/T (default: -0.1, cioè si aumenta la temparatura)
# echo "predictor" > $KSFN
DEBUG="1"
NP="1000"
alias awk='LANG=C awk'
EXES="../sim1statepnt_HC_MCNPTK.sh"
CNFCURI="startIso.cnf"
CNFCURN="startNem.cnf"
INIFILEI="startIso.cnf"
INIFILEN="startNem.cnf"
SCNF="start.cnf"
BFN="betaCur.dat"
BEFN="betaEnd.dat"
PFN="Pcur.dat"
KFN="EneVolKofke.dat"
KFNPRED="EneVolKofkePred.dat"
cp $KFN Ini-$KFN # backup
KSFN="kstat.dat"
NPI=`cat $INIFILEI| awk -F : '{if ($1=="parnum") print $2}'`
NPN=`cat $INIFILEN| awk -F : '{if ($1=="parnum") print $2}'`
EQSTEPS="50000"
PARF="ellipsoid_flex_mc.par"
#VINI_ISO=`cat VolIniIso.dat`
#VINI_NEM=`cat VolIniNem.dat`
BETAEND=`cat $BEFN`
#faccio così solo per fargli passare il while
#poi infatti betaCur.dat viene settato opportunamente
#in base al contenuto del file $KFN
LAST=`tail -1 $KFN`
BETAINI=`echo $LAST | awk '{print $1}'`
if [ "$1" == "" ]
then
DELBMOD="0.1"
DELB=`echo $BETAINI $BETAEND | awk -v delb=$DELBMOD '{if ($1 > $2) {printf("-%s",delb);} else {printf("%s",delb);}}'`
else
DELB="$1"
fi
fine="0"
echo "NPI= " $NPI " NPN=" $NPN " BETAINI=" $BETAINI " BETAEND=" $BETAEND " DELB=" $DELB
#exit 
while [ "$fine" == "0" ]
do 
LAST=`tail -1 $KFN`
PREV=`tail -2 $KFN | head -1`
PPREV=`tail -3 $KFN | head -1`
P0=`echo $LAST | awk '{print $2}'`
V0I=`echo $LAST | awk '{print $3}'`
U0I=`echo $LAST | awk '{print $4}'`
V0N=`echo $LAST | awk '{print $5}'`
U0N=`echo $LAST | awk '{print $6}'`
BETA0=`echo $LAST | awk '{print $1}'`
F0=`echo "-((${U0N})-(${U0I}) + (${P0})*(${V0N}-${V0I}))/(${BETA0}*(${P0})*(${V0N}-${V0I}))"|bc -l`
//echo "-((${U0N})-(${U0I}) + (${P0})*(${V0N}-${V0I}))/(${BETA0}*(${P0})*(${V0N}-${V0I}))"
//echo "F0="$F0
//exit
Pm1=`echo $PREV | awk '{print $2}'`
Vm1I=`echo $PREV | awk '{print $3}'`
Um1I=`echo $PREV | awk '{print $4}'`
Vm1N=`echo $PREV | awk '{print $5}'`
Um1N=`echo $PREV | awk '{print $6}'`
BETAm1=`echo $PREV | awk '{print $1}'`
Fm1=`echo "-((${Um1N})-(${Um1I}) + ${Pm1}*(${Vm1N}-${Vm1I}))/(${BETAm1}*${Pm1}*(${Vm1N}-${Vm1I}))"|bc -l`
Pm2=`echo $PREV | awk '{print $2}'`
Vm2I=`echo $PREV | awk '{print $3}'`
Um2I=`echo $PREV | awk '{print $4}'`
Vm2N=`echo $PREV | awk '{print $5}'`
Um2N=`echo $PREV | awk '{print $6}'`
BETAm2=`echo $PREV | awk '{print $1}'`
Fm2=`echo "-((${Um2N})-(${Um2I}) + ${Pm2}*(${Vm2N}-${Vm2I}))/(${BETAm2}*${Pm2}*(${Vm2N}-${Vm2I}))"|bc -l`
if [ "$KSTAT" == "predictor" ]
then
if [ "$PREV" == "" ]
then
#use trapezoid formula at begin
NEWP=`echo "${P0}*e((${DELB})*(${F0})*0.5)" | bc -l`
echo $NEWP > $PFN
NEWBETA=`echo "(${BETA0})+(${DELB})"| bc -l`
echo "$NEWBETA" > $BFN
else
#use midpoint formula 
NEWP=`echo "${Pm1}*e((${DELB})*(${Fm1})*2.0)" | bc -l`
echo $NEWP > $PFN
NEWBETA=`echo "${BETA0}+(${DELB})"| bc -l`
echo "$NEWBETA" > $BFN
fi
else
#CORRECTOR HERE
PPREV="`tail -3 $KFN | head -1`"
BETAp1=`echo $KFNPRED | awk '{print $1}'`
Pp1=`echo $KFNPRED | awk '{print $2}'`
Vp1I=`echo $KFNPRED | awk '{print $3}'`
Up1I=`echo $KFNPRED | awk '{print $4}'`
Vp1N=`echo $KFNPRED | awk '{print $5}'`
Up1N=`echo $KFNPRED | awk '{print $6}'`
Fp1=`echo "-((${Up1N})-(${Up1I}) + ${Pp1}*(${Vp1N}-${Vp1I}))/(${BETAp1}*${Pp1}*(${Vp1N}-${Vp1I}))"|bc -l`
if [ "$PREV" == "" ]
then
#use trapezoid formula at begin
NEWP=`echo "$P0*e(($DELB)*(($Fp1)+($F0))*0.5)" | bc -l`
echo $NEWP > $PFN
#NOTE: beta is the same as in predictor stage 
else
#use midpoint formula 
NEWP=`echo "${Pm1}*e(${DELB}*((${Fp1})+4.0*(${F0})+(${Fm1}))/3.0)" | bc -l`
echo $NEWP > $PFN
#NOTE: beta is the same as in predictor stage 
fi
fi
BETACUR=`cat $BFN`
ISOD="ISO-Beta-$BETACUR"
NEMD="NEM-Beta-$BETACUR"
EXEI="kofkeI-$BETACUR"
EXEN="kofkeN-$BETACUR"
RUNI=`ps ax | grep $EXEI`
RUNN=`ps ax | grep $EXEN`
PCUR=`cat $PFN`
KSTAT=`cat $KSFN`
if [ ! -e $ISOD ]
then
mkdir $ISOD
fi
cd $ISOD
if [ "$RUNI" == "" ] 
then
if [ \( -e "COORD_TMP0" \) -o \( -e "COORD_TMP1" \) ]
then
$EXEI -c >> screen
else
cp ../$CNFCURI $SCNF
TCUR=`echo "1.0/$BETACUR"| bc -l`
$EXES $PCUR $TCUR $EQSTEPS
fi
fi
cd ..
if [ ! -e $NEMD ]
then
mkdir $NEMD
fi
cd $NEMD
if [ "$RUNM" == "" ] 
then
if [ \( -e "COORD_TMP0" \) -o \( -e "COORD_TMP1" \) ]
then
$EXEM -c >> screen
else
cp ../$CNFCURN $SCNF
rm CorFinal
$EXES $PCUR $TCUR $EQSTEPS
fi
fi
cd ..
#wait for current iso and nem sims to finish here
if [ ! -e "$ISOD/CorFinal" ]
then
IS="running"
else
IS="finished"
fi
if [ ! -e "$NEMD/CorFinal" ]
then
NS="running"
else
NS="finished"
fi
while [ \( "$IS" == "running" \) -o \( "$NS" == "running" \) ]
do
if [ ! -e "$ISOD/CorFinal" ]
then
IS="running"
else
IS="finished"
fi
if [ ! -e "$NEMD/CorFinal" ]
then
NS="running"
else
NS="finished"
fi
sleep 30
done
#calculate running averages here and update Kofke file
N1="wc -l $ISOD/energy.dat"
FACT="echo "2.0/3.0"| bc -l"
AVENEISO=`cat $ISOD/energy.dat|awk -v NN=$NPI -v fact=$FACT -v n1=$N1 '{if (NR > fact*n1) {cc++; ene+=$2};} END { printf("%.15G\n",ene/cc/NN)}'`
N1="wc -l $NEMD/energy.dat"
AVENENEM=`cat $NEMD/energy.dat |awk -v NN=$NPN -v fact=$FACT -v n1=$N1 '{if (NR > fact*n1) {cc++; ene+=$2}; } END { printf("%.15G\n",ene/cc/NN)}'`
N1="wc -l $ISOD/volume.dat"
AVVOLISO=`cat $ISOD/volume.dat|awk -v NN=$NPI -v fact=$FACT -v n1=$N1 '{if (NR > fact*n1) {cc++; vol+=$2}; } END { printf("%.15G\n",vol/cc/NN)}'`
N1="wc -l $NEMD/volume.dat"
AVVOLNEM=`cat $NEMD/volume.dat|awk -v NN=$NPN -v fact=$FACT -v n1=$N1 '{if (NR > fact*n1) {cc++; vol+=$2}; } END { printf("%.15G\n",vol/cc/NN)}'`
#calculate new pressure here according to some integration scheme (e.g. trapezoid, midpoint, etc.)
#KFN:
#<beta> <pressure> <av. volume ISO> <av. energy ISO> <av. vol. NEM> <av. energy NEM> <#particles> 
if [ "$KSTAT" == "predictor" ]
then
echo "$BETACUR $PCUR $AVVOLISO $AVENEISO $AVVOLNEM $AVENENEM" > $KFNPRED
echo "corrector" > $KSTAT
else
echo "$BETACUR $PCUR $AVVOLISO $AVENEISO $AVVOLNEM $AVENENEM" >> $KFN
echo "predictor" > $KSTAT
fi
fine=`echo $BETAINI $BETAEND $BETACUR| awk '{if ($1 < $2) {if ($3 <= $2) printf("0"); else printf("1")}  else if ($3 <= $2) printf("1"); else printf("0");}'`
done
