#./create_GDNA_conf start_D2.cnf 2 13 28 0.45 2.0
#./create_GDNA_conf start_D3.cnf 2 10 21 0.45 3.0
# crea delle configurazioni con volume fraction uniformemente distribuite
# nell'intervallo 0.10-0.45 e le associa a pressioni ragionevolmente compatibili
# (stando alle equazioni di stato del Soft Matter 11, 2934 (2015))
NP="40"
BETA=`echo "1.0/0.124"| bc -l`
V0="113.097"
phi_1="0.1" # 0.1
phi_2="0.7"  #0.6
dphi=`echo "($phi_2-$phi_1)/$NP"| bc -l`
Pini="0.5"
DP="0.25" #"0.3"
cc="0"
while [ $cc -lt $NP ]
do
press=`echo $Pini $DP $cc $V0 $BETA | LANG=C gawk '{print ($1+$2*$3)/$4/$5}'`
phi=`echo $phi_1 $dphi $cc | LANG=C gawk '{print $1+$2*$3}'`
DIR="P_$press"
#FND2="start_D2_P${press}.cnf"
FND3="start_D3_P_${press}_cnf"
echo "PHI=" $phi "phi_1 =" $phi_1 " dphi=" $dphi
./create_ONEPATCH_conf $FND3 2 10 21 $phi 3.0 
cc=$[$cc+1]
done
