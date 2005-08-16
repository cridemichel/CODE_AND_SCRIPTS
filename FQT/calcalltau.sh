PERC=$HOME/ELLIPSOIDS/FQT/
if [ "$1" == "" ]
then
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
else
EL=$1
fi
FNT="alltautra-X0_${EL}.dat"
FNR="alltaurot-X0_${EL}.dat"
echo -n "" > $FNT
echo -n "" > $FNR
for f in Phi*
do 
cd $f
if [ ! -e ALL_DONE ]
then
cd ..
continue
fi
PHI=`echo $f | awk -F Phi '{print $2}'`
echo "Processing Phi=" $PHI
MAXQ=`$PERC/findmax sq0000.dat`
#echo "MAXQ TRA=" $MAXQ
TAU=`cat N-sqt.0000.k=$MAXQ | gawk 'BEGIN {xo=0.0; yo=0.0; K=1.0/exp(1.0);} {if ($2 < K) {m=($2-yo)/($1-xo); q=yo-m*xo; tau=(K-q)/m; printf("%f\n",tau); exit;} xo=$1; yo=$2; }'`
echo $EL $PHI $TAU >> $FNT
#echo "TAUT=" $TAU
MAXQ=`$PERC/findmax sq2020.dat`
#echo "MAXQ ROT=" $MAXQ
TAU=`cat N-sqt.2020.k=$MAXQ | gawk 'BEGIN {xo=0.0; yo=0.0; K=1.0/exp(1.0);} {  if ($2 < K) {m=($2-yo)/($1-xo); q=yo-m*xo; tau=(K-q)/m; printf("%f\n",tau); exit;} xo=$1; yo=$2; }'`
echo $EL $PHI $TAU >> $FNR
#echo "TAUR=" $TAU
cd ..
done
