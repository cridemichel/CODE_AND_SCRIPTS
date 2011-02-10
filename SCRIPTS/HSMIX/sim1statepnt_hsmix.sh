# $1 = volume fraction (Phi1)
# $2 = cicli di produzione (1 ciclo = 1 tempo di equilibratura) (se < 0 non fa la produzione)
# $3 = se < 0 non fa l'equilibratura, se = 0 usa il MSD per terminare l'equilibratura, se > 0 fa 
# $3 passi di equilibratura ( 0 = default ) 
# $4 = Phi2 o niente ed allora usa il nome della directory
# $5 = 1 cresce in due stadi con vol. frac. intermedia
#MR="mosrun"
MR=""
if [ "$1" = "" ]
then
echo "Syntax: sim1statepnt <Volume_Fraction Phi2> <production_cycles> [equilibration_steps]"
exit 0
fi
if [ "$2" = "" ]
then
echo "You have to supply the number of production cycles"
exit 0
fi
if [ "$3" = "" ]
then
EQSTPS=0
else
EQSTPS=$3
fi
PARFILE="silica_O1.par"
PHI1="$1"
if [ "$4" == "" ]
then
PHI2=`pwd | awk -F _ '{print $2}'` 
else
PHI2="$4"
fi
PHITOT=`echo "${PHI1}+${PHI2}"| bc -l`
if [ "$4" != "" ]
then
if [ ! -e Phi2_$4 ]
then
mkdir Phi2_$4
fi
cp ../$PARFILE Phi2_$4
cd Phi2_$4
else
if [ ! -e Phi1_$1 ]
then
mkdir Phi1_$1
fi
cp ../$PARFILE Phi1_$1
cd Phi1_$1
fi
rm -f COORD_TMP*
rm -f Store-*
HSMEXE="../../sticky_O1"
SIMRA="hsmixRA_Phi1_${PHI1}_Phi2_${PHI2}"
SIMGR="hsmixGR_Phi1_${PHI1}_Phi2_${PHI2}"
SIMEQ="hsmixEQ_Phi1_${PHI1}_Phi2_${PHI2}"
SIMPR="hsmixPR_Phi1_${PHI1}_Phi2_${PHI2}"
INIFILE="start.cnf"
GENCONF="../../genconfHS_Si"
STORERATE="0.01"
PARNUMA=250
SIGA="1.0"
SIGB="0.2"
DT="0.05"
#RANDOMIZZAZIONE INIZIALE
cp $PARFILE rand_$PARFILE
#echo "L:" $INIL >> rand_$PARFILE
PI="3.14159265358979"
TARGL=`echo "e((1.0/3.0)*l(${SIGA}^3*$PARNUMA*$PI/6.0/${PHI1}))" | bc -l`
PARNUMB=`echo "6.0*${PHI2}*($TARGL/$SIGB)^3/${PI}" | bc -l | awk '{printf("%d",$1)}'`
echo "Phi1=" $PHI1 " PHI2=" $PHI2 " TARGL=" $TARGL " PARNUMB=" $PARNUMB
#exit
#genera la configurazione iniziale
echo "PARNUMA=" $PARNUMA " PARNUMB=" $PARNUMB
$GENCONF $PARNUMA $PARNUMB > $INIFILE
#>>> SET TEMPERATURE TO 1.0
../../set_params.sh rand_$PARFILE inifile $INIFILE stepnum 100 targetPhi 0.0 storerate 0.0 intervalSum 0.2 endfile ${SIMRA}.cor  
ln -sf $HSMEXE $SIMRA
./$SIMRA -fa ./rand_${PARFILE} > screen_$SIMRA 
#exit 
#echo "FINE SET TEMP"
#>>> EQUILIBRATE STARTING DENSITY (I.E. RANDOMIZE)
PARNUM=`echo "${PARNUMA}+${PARNUMB}" | bc`
MAXCOLL=`echo "${PARNUM}*200"| bc`
../../set_params.sh rand_$PARFILE stepnum 50000 maxcoll $MAXCOLL targetPhi 0.0 storerate 0.0 intervalSum 2.0 scalevel 0 inifile ${SIMRA}.cor endfile ${SIMRA}.cor 
ln -sf $HSMEXE $SIMRA 
$MR ./$SIMRA -f ./rand_${PARFILE} >> screen_$SIMRA 
#echo "FINE RANDOMIZZAZIONE"
#CRESCITA
if [ "$5" != "" ]
then
PHITOTINT=`echo "$PHITOT 0.9" | awk '{print ($1*$2)}'`
../../set_params.sh $PARFILE stepnum 10000000 maxcoll -1 targetPhi $PHITOTINT storerate 0.0 intervalSum 0.05 inifile ${SIMRA}.cor endfile ${SIMGR}.cor 
ln -sf $HSMEXE $SIMGR
$MR ./$SIMGR -f ./$PARFILE > screen_${SIMGR}_int
MAXCOLL=`echo "${PARNUM}*1000"| bc`
../../set_params.sh $PARFILE stepnum 1000000 maxcoll $MAXCOLL targetPhi 0.0 storerate 0.0 intervalSum 2.0 inifile ${SIMGR}.cor endfile ${SIMRA}.cor 
ln -sf $HSMEXE $SIMGR
$MR ./$SIMGR -f ./$PARFILE > screen_${SIMGR}_inteq 
#PHITOTINT=`echo "$PHITOT 0.95" | awk '{print ($1*$2)}'`
#../../set_params.sh $PARFILE stepnum 10000000 maxcoll -1 targetPhi $PHITOTINT storerate 0.0 intervalSum 0.01 inifile ${SIMRA}.cor endfile ${SIMGR}.cor 
#ln -sf $HSMEXE $SIMGR
#$MR ./$SIMGR -f ./$PARFILE > screen_${SIMGR}_int
#MAXCOLL=`echo "${PARNUM}*1000"| bc`
#../../set_params.sh $PARFILE stepnum 1000000 maxcoll $MAXCOLL targetPhi 0.0 storerate 0.0 intervalSum 2.0 inifile ${SIMGR}.cor endfile ${SIMRA}.cor 
ln -sf $HSMEXE $SIMGR
#$MR ./$SIMGR -f ./$PARFILE > screen_${SIMGR}_inteq 
fi
../../set_params.sh $PARFILE stepnum 10000000 maxcoll -1 targetPhi $PHITOT storerate 0.0 intervalSum 0.01 inifile ${SIMRA}.cor endfile ${SIMGR}.cor 
ln -sf $HSMEXE $SIMGR
$MR ./$SIMGR -f ./$PARFILE > screen_$SIMGR 
#exit 
#EQUILIBRATURA
if [ $EQSTPS -ge 0 ]
then
if [ $EQSTPS -eq 0 ]
then
STPS=10000000
#soglia messa a 2 x sigma
TMSDA=2.0
TMSDB=0.4
else
STPS=$EQSTPS
TMSDA="-1.0"
TMSDB="-1.0"
fi
../../set_params.sh $PARFILE stepnum $STPS targetPhi 0.0  storerate 0.0 intervalSum 5.0 tmsd2endA $TMSDA tmsd2endB $TMSDB inifile ${SIMGR}.cor endfile ${SIMEQ}.cor
ln -sf $HSMEXE $SIMEQ 
$MR ./$SIMEQ -f ./$PARFILE > screen_$SIMEQ 
fi
#PRODUZIONE
if [ $2 -gt 0 ]
then
if [ $[EQSTPS] -gt 0 ]
then
STCI="$EQSTPS"
STPS=`echo "$STCI*$2"| bc -l`
#echo "qui STPS=$STPS"
else
STCI=`cat screen_$SIMEQ | awk '{if ($1=="[MSDcheck]") print $3}'`
STPS=`echo "$STCI*$2"| bc -l`
fi
NN=`echo "l($DT*$STCI/$STORERATE)/l(1.3)" | bc -l | awk '{printf("%d",$0)}'`
../../set_params.sh $PARFILE stepnum $STPS targetPhi 0.0 storerate $STORERATE intevalSum 5.0 tmsd2endA -1.0 tmsd2endB -1.0 NN $NN inifile ${SIMEQ}.cor endfile ${SIMPR}.cor
ln -sf $HSMEXE $SIMPR
$MR ./$SIMPR -f ./$PARFILE > screen_$SIMPR 
fi
cd ..
