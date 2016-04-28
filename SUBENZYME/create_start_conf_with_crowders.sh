#UI=("5.0" "5.0" "10.0" "10.0")
#UO=("10.0" "5.0" "5.0" "10.0")
UI=("10.0")
UO=("10.0")
SIGS="1.0"
SIGE="4.5"
SIGC="10.2"
NENZ="338" # formerly NENZ=338
#PPL="45"
PPL="61"
NCROWARR=("160000" "133000" "107000" "80000" "53000" "27000")
NCROW=${NCROWARR[0]} 
#NSUB=`echo ${PPL}*${PPL}*${PPL}-$NENZ|bc`
NSUB=`echo ${PPL}*${PPL}*${PPL}-$NENZ-$NCROW|bc`
#NSUB=`echo 45*45*45-$NENZ|bc`
i=0
while [ $i -lt 1 ]
do
#startCROW_phi0.40.cnf 15 15 15 0.40 4.5 1.0 9.4 4.5 1.0 9.4 1 1000
#PHI="0.20"
LBOX="670.080659652254"
echo "UI=" ${UI[$i]} " UO=" ${UO[$i]}
./create_SUBENZ_conf start_CROWD_uI_${UI[$i]}_uO_${UO[$i]}.cnf $PPL $PPL $PPL ${LBOX} $SIGE $SIGS $SIGC $SIGE $SIGS $SIGC $NENZ $NSUB ${UI[$i]} ${UO[$i]}
i=$[i+1]
done
