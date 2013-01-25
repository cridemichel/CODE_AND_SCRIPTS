export LC_NUMERIC=C
FN="all_pts.dat"
DELOSF="iggres_vs_sig.dat"
RA="eqconst.dat"
#numero di Avogadro
Nav="6.023E23"
#conversione da 1 unità ridotta di lunghezza a 1 m
ru2m="2.0E-09"
#il secondo argomento è la frazione di punti da scartare, ad es. 0.6 vuol dire che scarterà il 60% dei punti partendo dal primo
if [ "$2" = "" ]
then
DN="half"
fi
echo -n "" > $FN
echo -n "" > $DELOSF
for f in `cat $1`
do
cd $f
echo "DIR=" $f
NP=`wc -l bi-mono-bonds.dat | awk '{print $1}'`
if [ "$DN" = "half" ]
then
NP=`echo "$NP/2" | bc `
else
NP=`echo "$NP*$2"|bc -l|awk '{printf ("%d",$1)}'` 
fi
echo "Points to skip #" $NP
L=`tail -1 CorIni`
NANT=`cat CorIni | awk '{if (nat==1) {print $5; nat=2}; if ($0=="@@@") nat++;}'`
V=`echo "$L*$L*$L" | bc -l`
S=`echo "$L*$L" | bc -l`
CF=`echo "$V*1000*$ru2m^3*$Nav"|octave -q| awk -F '=' '{print $2}'`
CFS=`echo "$S*$ru2m^2*$Nav"|octave -q| awk -F '=' '{print $2}'`
BINI=`echo $f | awk -F '/' '{print $1}' | awk -F '_' '{print $2}'`
SIG=`echo $f | awk -F '/' '{print $2}' | awk -F '_' '{print $2}'`
#echo "CF=" $CF "CFS=" $CFS
echo "NP="$NP "CF=" $CF " CFS=" $CFS "NANT=" $NANT
echo $BINI $SIG `cat bi-mono-bonds.dat | LC_NUMERIC=C awk -v np="$NP" -v nant="$NANT" -v cf="$CF" -v cfs="$CFS" '{if (NR >= np) {B+=$2; sigB+=$2*$2; mono+=$3; bi+=$4; cc++}} END{print (B/cc/cf, mono/cc/cfs, bi/cc/cfs, (nant-(mono/cc)-(bi/cc)*2)/cfs, (bi/cc+mono/cc)/nant, nant/cfs, sqrt(sigB/cc-B*B/(cc*cc))/(B/cc),mono/cc/250,bi/cc/250)}'` >> ../../$FN
echo $SIG `cat bi-mono-bonds.dat | LC_NUMERIC=C awk -v np="$NP" -v nant="$NANT" -v vol="$V" -v sur="$S" 'BEGIN {BINI=-1} {if (NR >= np) {if (BINI==-1) BINI=$2+$3+$4;  B+=$2; sigB+=$2*$2; mono+=$3; bi+=$4; cc++}} END{print (BINI, B/cc, mono/cc+bi/cc, mono/cc, bi/cc, BINI/vol, BINI/sur, (mono/cc+bi/cc)/sur)}'` >> ../../$DELOSF
cd ..
cd ..
done
cat $FN | sort -k 1 -n > sorted$FN
cat $DELOSF|sort -k 1 -n > sorted-$DELOSF
#cat sorted$FN | LC_NUMERIC=C awk -v LA="$L" -v cf="$CF" -v cfs="$CFS" -v sup="$S" '{D=1-2*33.5*$3*cf/(LA*LA*LA); K1=$4/($3*$6)/(1.0-2.0*33.5*$3*cf/(LA*LA*LA))/(1-0.0*(30.0*$4*cfs/sup-64.0*$5*cfs/sup)); K2=$5/($4*$6*(1-exp(-3.1415*144*$2))); print ($8,K1,K2,$8,$3,D)}' > $RA
cat sorted$FN | LC_NUMERIC=C awk -v LA="$L" -v cf="$CF" -v cfs="$CFS" -v sup="$S" '{D=1-2*33.5*$3*cf/(LA*LA*LA); K1=$4/($3*$6); K2=$5/($4*$6*(1-exp(-3.1415*144*$2))); print ($8,K1,K2,$8,$3,D)}' > $RA
