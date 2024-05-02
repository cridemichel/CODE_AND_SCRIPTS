TO_ADD="DQtensSteps: 100\n\
DQtensCalc: 100\n\
DQtensName: DQtens-\n\
PtensSteps: 100\n\
PtensCalc: 100\n\
PtensName: Ptens-"
for d in `ls -d RUN_*`
do
OD=`pwd`
cd $d
 echo $TO_ADD >> ellipsoid_dGEBA.par 
cd $OD
done
