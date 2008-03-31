PFQT=$HOME/ELLIPSOIDS/FQT/
PMSD=$HOME/ELLIPSOIDS/MSD/
QB=5
QE=30
PTS=100000
LOGFILE=analysis.log
if [ "$1" == "" ] 
then
TEMPS=`ls -d T*`
else
TEMPS="$1"
fi
for f in $TEMPS
do
cd $f
echo -n "Dir:" $f "..."
echo -n "" > $LOGFILE 
rm -fr RHOTMPA/ RHOTMPB/
ls Cnf* | sort -t - -k 3 -n > lista
cat lista | awk '{if (NR % 20 == 1) print $0}' > listaSq
$PFQT/calcSq -c listaSq $QB $QE
QMAX=`$PFQT/findmax Sq.dat 0`
$PFQT/calcrho lista $QMAX $QMAX >> $LOGFILE 2>&1
$PFQT/calcfqtcoll --ncomps 2 $QMAX $QMAX $PTS >> $LOGFILE 2>&1 
$PFQT/calcfqtself lista $PTS $QMAX $QMAX >> $LOGFILE 2>&1
$PMSD/calcmsd -ng lista $PTS >> $LOGFILE 2>&1
cd ..
echo "done!"
done
