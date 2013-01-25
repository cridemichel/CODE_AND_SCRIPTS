ROOTDIR=".."
ls -d ISO-Beta-* | sort -t - -k 3 -n > lista
for f in `cat lista`
do
BBEG=`echo $f| awk -F - '{print $3}'`
BFIN=`echo $BBEG+0.3| bc -l| awk '{printf("%.1f\n",$0)}'`
ISOINI="ISO-Beta-$BBEG"
NEMINI="NEM-Beta-$BBEG"
if [ \( ! -e $ISOINI/CorFinal \) -o \( ! -e $NEMINI/CorFinal \) ]
then
continue
fi
DN="BETA-${BBEG}-${BFIN}-ELONG-h-0.1-NEW"
BETA=`echo $f | awk -F - '{print $3}'`
if [ ! -e $ROOTDIR/$DN ]
then
mkdir $ROOTDIR/$DN
fi
cp $ISOINI/CorFinal $ROOTDIR/$DN/startIso.cnf
cp $NEMINI/CorFinal $ROOTDIR/$DN/startNem.cnf
cp ellipsoid_flex_mc.par $ROOTDIR/$DN
cp set_*                 $ROOTDIR/$DN
cp ../KOFKE_h0.1_SCRIPTS/kofke_run.sh $ROOTDIR/$DN
cp init_kofke.sh         $ROOTDIR/$DN
cp sim1statepnt_HC_MCNPTK.sh $ROOTDIR/$DN
cp EneVolKofke.dat       $ROOTDIR/$DN
cp ellipsHC              $ROOTDIR/$DN
CD=`pwd`
cd $ROOTDIR/$DN
cat EneVolKofke.dat | awk -v bbeg="$BBEG" '{if ($1==bbeg) print $0}' > _a_a_a_
cp _a_a_a_ EneVolKofke.dat
cp EneVolKofke.dat EneVolKofke_bak.dat
echo "predictor" > kstat.dat
echo $BFIN > betaEnd.dat
rm _a_a_a_
cd $CD
done
