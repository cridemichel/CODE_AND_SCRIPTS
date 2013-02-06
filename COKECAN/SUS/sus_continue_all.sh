for f in N_*
do 
cd $f
# FIND THE LATEST COORD_TMP_ASCII FILE (and check that it is not corrupted)
CFO=`ls -rt COORD_TMP_ASCII*| head -1`
CFN=`ls -rt COORD_TMP_ASCII*| tail -1`
PNO=`cat $CFO| awk -F : '{if ($1=="parnum") {print $2; exit;}}'`
PNN=`cat $CFN| awk -F : '{if ($1=="parnum") {print $2; exit;}}'`
POSATO=`cat $CFO | awk 'BEGIN{at=0} {if ($0=="@@@") {at++;}; if (at==3) {print NR; exit;}}'`
POSATN=`cat $CFN | awk 'BEGIN{at=0} {if ($0=="@@@") {at++;}; if (at==3) {print NR; exit;}}'`
WCO=`cat $CFO| wc -l`
WCN=`cat $CFN| wc -l`
EXPLO=`echo $POSATO+$PNO+1|bc`
EXPLN=`echo $POSATN+$PNN+1|bc`
echo "WCO=" $WCO " CFO=" $CFO " PNO=" $PNO " POSATO=" $POSATO " EXPLO=" $EXPLO
echo "WCN=" $WCN " CFN=" $CFN " PNN=" $PNN " POSATN=" $POSATN " EXPLN=" $EXPLN
if [ $[EXPLN] -eq $[WCN] ]
then
echo "I chose CFN=" $CFN
CF="$CFN" 
elif [ $[EXPLO] -eq $[WCO] ]
then
echo "I chose CFO=" $CFO
CF="$CFO" 
else
echo "DIR=" $DN
echo "FATAL ERROR: both COORD_TMP_ASCII corrupted!"
echo "I skip this directory..."
cd ..
continue
fi
nohup mosrun ./HC_$f -ca $CF > screen &
sleep 1
cd ..
done
