if [ "$3" != "" ]
then
MINMAX="1"
else
MINMAX="0"
fi
FNT="rates_vs_S.dat"
echo -n "" > $FNT
if [ "$1" == "" ]
then 
echo "expfit.sh <lista_dirs>"
exit
fi
if [ "$MINMAX" == "1" ]
then
TMIN="$2"
TMAX="$3"
SMIN=`cat $1 | LANG=C gawk -F _ 'BEGIN {MIN=-1;}{if (MIN==-1 || $2 < MIN) MIN=$2;} END {print MIN}'`
SMAX=`cat $1 | LANG=C gawk -F _ 'BEGIN {MAX=-1;}{if (MAX==-1 || $2 > MAX) MAX=$2;} END {print MAX}'`
echo "SMIN= " $SMIN " SMAX= " $SMAX
echo "TMIN= " $TMIN " TMAX= " $TMAX
TMMAX="$4"
fi
for f in `cat $1` 
do
cd $f
#echo "f=" $f " pwd=" `pwd`
CMAX=`cat startGR.cnf| LANG=C gawk 'BEGIN {at=0}; {if ($0=="@@@" && at==0) {at=1; nra=NR+1;} if (at==1 && NR==nra) print $1;}'` # awk '{if (NR==N+2) print $2; if ($1=="@@@") N=NR;}'`
#CMAX=`cat startEQ.cnf| LANG=C gawk '{if (NR==7) print $1}'` 
CMIN="0"
echo "CMAX= " $CMAX
MS=`echo $f| awk -F _ '{print $2}'`
if [ "$MINMAX" == "1" ]
then
STA=`echo $TMAX $TMIN $SMAX $SMIN| LANG=C gawk -v s=$MS -v tmin=$TMIN -v tmax=$TMAX -v smax=$SMAX '{print tmin-(tmin-tmax)*(s-$3)/($4-$3)}'`
NRI=`wc -l SP-popul.dat | LANG=C gawk '{printf("%d",2*$1/3)}'`
CAVG=`cat SP-popul.dat | LANG=C gawk -v nrI=$NRI '{if (NR > nrI) {sum+=$4; NP+=1;}} END {print (sum/NP)}'`
STA=`echo $TMAX $TMIN $CMAX $CMIN| LANG=C gawk -v s=$CAVG -v tmin=$TMIN -v tmax=$TMAX '{print tmin-(tmin-tmax)*(s-$3)/($4-$3)}'`
STB=`echo $STA $TMMAX | awk '{print $1*$2}'`
#STB="$TMMAX"
fi
echo "CAVG= " $CAVG
#echo "STA= " $STA " STB=" $STB
#exit
echo "Processing file=" $f
FN="SP-popul.dat"
#considera tutti valori tra STAF e STBF del tempo a cui è arrivato
if [ "$MINMAX" == "0" ]
then
STAF="0.05"
STBF="0.20"
TF=`tail -1 $FN | awk '{print $1}'`
STA=`echo ${STAF} $TF | LANG=C gawk '{print $1*$2}'` #0.001
STA="1000"
fi
L=`tail -n 1 startGR.cnf| awk '{print $1}'`
PN=`cat startGR.cnf | grep parnum | LANG=C gawk -F : '{print $2}'`
echo "L=" $L " PN=" $PN
if [ "$MINMAX" == "0" ]
then
STB=`echo ${STBF} $TF | LANG=C gawk '{print $1*$2}'` #`cat $FN| tail -1 | LANG=C awk '{print 9.*$1/10.}'`
STB="50000"
fi
echo "===> STA=" $STA " STB=" $STB 
A1="100000.0"
B1=100.0
C1=100.0
NS=`cat startGR.cnf| LANG=C gawk 'BEGIN {at=0}; {if ($0=="@@@" && at==0) {at=1; nra=NR+1;} if (at==1 && NR==nra) print $2;}'` # awk '{if (NR==N+2) print $2; if ($1=="@@@") N=NR;}'`
#NS=`cat startGR.cnf| LANG=C gawk '{if (NR==7) print $2}'` 
NE=`cat startGR.cnf| LANG=C gawk 'BEGIN {at=0}; {if ($0=="@@@" && at==0) {at=1; nra=NR+1;} if (at==1 && NR==nra) print $1;}'` # awk '{if (NR==N+2) print $2; if ($1=="@@@") N=NR;}'`
#NE=`cat startGR.cnf| LANG=C gawk '{if (NR==7) print $1}'` 
S0=`echo $NS $L | LANG=C gawk '{print $1/$2/$2/$2}'`
E0=`echo $NE $L | LANG=C gawk '{print $1/$2/$2/$2}'`
DS=`echo "0.1 2.0"| LANG=C gawk '{printf("%f", $1/$2)}'`
RR=`echo "4.5 1.0 2.0" | LANG=C gawk '{printf("%f", ($1+$2)/$3)}'`
#echo "5000.0 $L"|awk '{print($1/($2*$2*$2))}*$2*$2'
#RHO=`echo "5000.0 ${L} ${L} ${L}"| LANG=C gawk '{printf("%G", ($1/($2*$2*$2)))}'`
#KKR=`echo ${B1} ${S0} ${E0} |LANG=C gawk '{print $2/$1/$3}'`
AVGS=`cat SP-popul.dat | LANG=C gawk -v tmin=$STA -v tmax=$STB '{if ($1 > tmin && $1 < tmax) {avgs+=$2; cc++;}} END\\
 {print avgs/cc}'`
AVGE=`cat CSE_vs_t.dat | LANG=C gawk -v tmin=$STA -v tmax=$STB -v ne=$NE '{if ($1 > tmin && $1 < tmax) {avge+=(ne-$2); cc++;}} END\\
 {print avge/cc}'`
AVGES=`cat CSE_vs_t.dat | LANG=C gawk -v tmin=$STA -v tmax=$STB -v ne=$NE '{if ($1 > tmin && $1 < tmax) {avge+=($2); cc++;}} END {print avge/cc}'`
AVGR=`cat ratesSE.dat | LANG=C gawk -v tmin=$STA -v tmax=$STB 'BEGIN {first=0;} {if (first==1){tini=$1; \\
nevini=$4*$1; first=2;} if ($1 > tmin && first==0) first=1; if ($1 > tmax) {print ($4*$1-nevini)/($1-tini); exit;} }'`
AVGKS=`cat ratesSE.dat | LANG=C gawk -v tmin=$STA -v tmax=$STB 'BEGIN {first=0;} {if (first==1){tini=$1; \\
nevini=$2*$1; first=2;} if ($1 > tmin && first==0) first=1; if ($1 > tmax) {print ($2*$1-nevini)/($1-tini); exit;} }'`
AVGKMD=`cat ratesSE.dat | LANG=C gawk -v tmin=$STA -v tmax=$STB 'BEGIN {first=0;} {if (first==1){tini=$1; \\
nevini=$3*$1; first=2;} if ($1 > tmin && first==0) first=1; if ($1 > tmax) {print ($3*$1-nevini)/($1-tini); exit;} }'`
echo "AVGE= " $AVGE " AVGR=" $AVGR
###KKR=`echo ${B1} ${S0} ${E0} |LANG=C gawk '{print $2/$1/$3}'`
#KKR=`echo ${B1} ${S0} ${E0} |LANG=C gawk '{print $1/$3}'`
#KKR=`echo ${B1} ${MS}|awk '{print $2/$1}'` # così uso la concentrazione molare di S (MS) e non S0 (number density)
KD=`echo $AVGR $AVGE $AVGS $L| LANG=C gawk '{print $4^6*$1/$2/$3}'`
KS=`echo $AVGKS $AVGES $L | LANG=C gawk '{print $3^3*$1/$2}'`
KMD=`echo $AVGKMD $AVGES $L | LANG=C gawk '{print $3^3*$1/$2}'`
echo "PP=" $PP " L=" $L " RR=" $RR  " DS=" $DS  " NS= " $NS " S0=" $S0  " E0=" $E0 " CMIN= "$CMIN " CMAX=" $CMAX
#KK=`echo "($PN/$L/$L/$L)/${B1}" | bc -l`
#KK=`echo "100000/${B1}/${SIG}/($PN/$L/$L/$L)" | bc -l`
#KK=`echo "100000/${B1}/(${SIG}+1.0)/(1.0+1.0/${SIG})/($PN/$L/$L/$L)" | bc -l`
#echo "$MS $KKR $KK $RHO $PP"  >> ../$FNT
echo "$S0 $KS $KMD $KD"  >> ../$FNT
###
#rm fit.tmp
cd ..
done
