rm -fr *-Beta-*
echo "predictor" > kstat.dat
cp EneVolKofke_bak.dat EneVolKofke.dat
P=`cat kofke_run_PID`
if [ "$P" != "" ]
then
kill $P
fi
rm kofke_run_PID
