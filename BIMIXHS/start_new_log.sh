# $1=temp $2=NN
PHI="0.09"
cp ../bimixhsSWBar.par .
cat bimixhsSWBar.par | awk -v temp=$1 -F ':' '{if ($1=="temperat") print ("temperat: ",temp); else print $0}' > _aaa_
cp _aaa_ bimixhsSWBar.par
cat bimixhsSWBar.par | awk -v NN=$2 -F ':' '{if ($1=="NN") print ("NN: ",NN); else print $0}' > _aaa_
cp _aaa_ bimixhsSWBar.par
rm _aaa_
CD=`pwd`
cd ../T$1
LASTCNF=`ls Cnf* | sort -t - -k 2 -n | tail -1`
cd $CD 
cp ../T$1/$LASTCNF .
cp ../T$1/fra2cri.sh .
./fra2cri.sh $LASTCNF > pwm.res
rm $LASTCNF
ln -sf $HOME/BIMIXHS/MDsimul/bin/bimixhs bmT$1PHI$PHI 
./bmT$1PHI$PHI -fa bimixhsSWBar.par > screen &
