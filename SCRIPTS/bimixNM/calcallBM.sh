PFQT=$HOME/ELLIPSOIDS/FQT
PMSD=$HOME/ELLIPSOIDS/MSD
QB=20
QE=20
PTS=1000
LOGFILE=analysis.log
for f in T*
do
cd $f
echo -n "Dir:" $f "..."
echo -n "" > $LOGFILE 
rm -fr RHOTMPA/ RHOTMPAB/
ls Cnf* | sort -t - -k 3 -n > lista
$PFQT/calcrho lista $QB $QE >> $LOGFILE 2>&1
$PFQT/calcfqtcoll $QB $QE $PTS >> $LOGFILE 2>&1 
$PFQT/calcfqtself lista $PTS $QB $QE >> $LOGFILE 2>&1
$PMSD/calcmsd lista $PTS >> $LOGFILE 2>&1
cd ..
echo "done!"
done
