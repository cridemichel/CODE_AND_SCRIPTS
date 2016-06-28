#./create_GDNA_conf start_D2.cnf 2 13 28 0.45 2.0
#./create_GDNA_conf start_D3.cnf 2 10 21 0.45 3.0
# crea delle configurazioni con volume fraction uniformemente distribuite
# nell'intervallo 0.10-0.45 e le associa a pressioni ragionevolmente compatibili
# (stando alle equazioni di stato del Soft Matter 11, 2934 (2015))
NP="40"
phi_1="0.1" # 0.1
phi_2="0.6"  #0.6
dphi=`echo "($phi_2-$phi_1)/$NP"| bc -l`
T="0.124"
DIAM="3.0"
LEN="16.0"
PATCHDIAM="4.0"
VOL=`echo $DIAM $LEN| LANG=C gawk '{print 3.14159*$2*$1/2.0}'`
echo "VOL=$VOL"
PFACT=`echo $T $VOL | LANG=C gawk '{print $2/$1}'`
Pini=`echo 0.5 $PFACT| LANG=C gawk '{print $1*$2}'`
DP=`echo 0.2 $PFACT | LANG=C gawk '{print $1*$2}'`
cc="0"
while [ $cc -lt $NP ]
do
press=`echo $Pini $DP $cc | LANG=C gawk '{print $1+$2*$3}'`
phi=`echo $phi_1 $dphi $cc | LANG=C gawk '{print $1+$2*$3}'`
DIR="P_$press"
#FND2="start_D2_P${press}.cnf"
FND3="start_D3_P${press}.cnf"
echo "PHI=" $phi "phi_1 =" $phi_1 " dphi=" $dphi
#./create_HCBIGSPOT_conf $FND2 2 13 28 $phi 2.0
#./create_HCBIGSPOT_conf $FND3 2 10 21 $phi 3.0 
./create_HCBIGSPOT_conf $FND3 4 20 22 $phi $DIAM $PATCHDIAM
cc=$[$cc+1]
done
