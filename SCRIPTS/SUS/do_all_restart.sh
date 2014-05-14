TYPE="restart" 
#TYPE="start"
if [ $2 = "" ]
then
echo "you must provide Nbeg and Nfin"
exit
fi
Nbeg=$1
Nfin=$2
if [ $3 == "" ]
then
BS="100"
else
BS="$3"
fi
echo $BASHPID > do_all_restart.sh.PID
rm -f N_*_*/screen
if [ $TYPE = "start" ]
then
./start_sus_para.sh $Nbeg $BS 40
else
./restart_sus_para.sh $Nbeg $BS 40
fi
sleep 120
N=$BS
F="no"
while [ $F = "no" ]
do
R=`./count_running.sh`
T=`./count_terminated.sh`
if [ $T -lt $Nfin ]
then
if [ $R -lt $BS ]
then
echo "N=" $N "/ terminated # " $T " / running # " $R
echo "Starting new simulation"
if [ $TYPE = "start" ]
then
./start_sus.sh $N $[$N+1]
else
./restart_sus.sh $N $[$N+1]
fi
N=$[$N+1]
fi
else
F="yes"
fi
if [ $N -ge $Nfin ]
then
F="yes"
fi
sleep 5
done
