rm -fr *-Beta-*
echo "predictor" > kstat.dat
cp EneVolKofke_bak.dat EneVolKofke.dat
P=`cat kofke_run_PID`
kill $P
rm kofke_run_PID
