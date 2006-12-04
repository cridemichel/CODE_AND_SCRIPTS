#!/bin/bash
#normalmente se i file ci sono già non ricalcola le varie grandezze ma si puo' forzare il ricalcolo
#alha stimato dal caso X0=0.33333333 Phi=0.45
#(ved. fig. tau_vs_q_X0_0.33333333_Phi0.45.agr)
ALPHA='1.6' 
EXE_PATH=$HOME/postdoc/hardellipsoid/hardellSVN/CODE/
FN="all_qmax.dat"
echo -n "" > $FN
for f in X0*
do 
echo "Doing " $f
cd $f
for g in Phi*
do
cd $g
if [ ! -e Sq.dat ]
then
echo $f " " $g " Sq.dat does not exist!"
cd ..
continue
fi
#getting L
if [ \( -e "Store-0-0.gz" \) -o \( -e "Store-0-0" \) ]
then
if [ -e "Store-0-0.gz" ]
then
gunzip Store-0.0.gz
RIZIP=1
else
RIZIP=0
fi
L=`tail -1 Store-0-0` 
if [ "$RIZIP" == "1" ]
then
gzip Store-0-0
fi
else
if [ -e "rand_ellips.par" ]
then
L=`cat rand_ellips.par | awk -F ':' '{ if ($1=="L") print $2}'` 
else
echo "I can not get L! Skipping state point" $f " " $g
cd ..
continue
fi
fi
echo "L= " $L " Doing " $g
PR=$EXE_PATH/FQT/findmax
SQMAX=`$PR Sq.dat 0`
#PI=`echo "pi" |octave| LANG=C gawk -F '=' '{if ($1=="ans") print $2}'`
SQMAXFQS=`echo "$SQMAX" | LANG=C gawk -v sqmax="$SQMAX" -v L="$L" -v alpha="$ALPHA" '{pi=atan2(1.0,1.0)*4.0; printf("%.15G",((2.0*pi/L)*(0.5*sqmax+1.25))^alpha)}'`
X0=`echo $f | gawk -F '_' '{print $2}'`
PHI=`echo $g | LANG=C gawk -F 'Phi' '{print $2}'`
echo $X0 " " $PHI " " $SQMAXFQS >> ../../$FN
cd ..
done
cd ..
done
