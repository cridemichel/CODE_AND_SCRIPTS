FN="all_pts.dat"
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
for f in $1
do
cd $f
NP=`wc -l bi-mono-bonds.dat | awk '{print $1}'`
if [ "$DN" = "half" ]
then
NP=`echo "$NP/2" | bc `
else
NP=`echo "$NP*$2"|bc -l|awk '{printf ("%d",$1)}'` 
fi
echo "Points to skip #" $NP
L=`tail -1 CorIni`
V=`echo "$L*$L*$L" | bc -l`
S=`echo "$L*$L" | bc -l`
CF=`echo "1000*$ru2m^3*$Nav"|octave -q| awk -F '=' '{print $2}'`
CFS=`echo "$ru2m^2*$Nav"|octave -q| awk -F '=' '{print $2}'`
BINI=`echo $f | awk -F '/' '{print $1}' | awk -F '_' '{print $2}'`
SIG=`echo $f | awk -F '/' '{print $2}' | awk -F '_' '{print $2}'`
echo "NP="$NP "CF=" $CF " CFS=" $CFS
echo $BINI $SIG `cat bi-mono-bonds.dat | LANG=C awk -v np="$NP" -v cf="$CF" -v cfs="$CFS" '{if (NR >= np) {B+=$2; mono+=$3; bi+=$4; cc++}} END{print (B/cc/cf, mono/cc/cfs, bi/cc/cfs)}'` >>../../$FN
cd ..
cd ..
done
