UI=("5.0" "5.0" "10.0" "10.0")
UO=("10.0" "5.0" "5.0" "10.0")
SIGS="1.0"
SIGE="4.5"
SIGC="9.4"
i=0
while [ $i -lt 4 ]
do
#startCROW_phi0.40.cnf 15 15 15 0.40 4.5 1.0 9.4 4.5 1.0 9.4 1 1000
PHI="0.0000536857"
echo "UI=" ${UI[$i]} " UO=" ${UO[$i]}
./create_SUBENZ_conf start_NOCROWD_uI_${UI[$i]}_uO_${UO[$i]}.cnf 30 30 30 ${PHI} $SIGE $SIGS $SIGC $SIGE $SIGS $SIGC 338 26662 ${UI[$i]} ${UO[$i]}
i=$[i+1]
done
