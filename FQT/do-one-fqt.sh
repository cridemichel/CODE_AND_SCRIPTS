cd Phi$1
if [ -e IGNORE_THIS ]
then 
echo "IGNORING AS YOU REQUEST"
cd ..
exit
fi
if [ -e ALL_DONE ]
then 
echo "ALL DONE HERE!"
cd ..
exit 
fi
touch IN_PROGRESS
gunzip Store*gz
ls Store* > junk
cat junk | while read nomefile
do
echo  $nomefile
awk -F":" 'BEGIN {cc=0} {if (cc>0) print $0; if ($1=="@@@") cc=cc+1; }' $nomefile > Cnf$nomefile
done
gzip Store*
PERC=$HOME/ELLIPSOIDS/FQT/
ls Cnf* | sort -t - -k 2 -k 3 -n > listaconf
q=2
while [ $q -lt 100 ]
do
$PERC/rho-dip << !
$q
!
q=$[$q+1]
done
rm -f CnfStore* 
PR=$HOME/ELLIPSOIDS/FQT/calcfqt
if [ "$1" == "" ]
then
NN=`cat Store-0-0.gz | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
else
NN=$1
fi
echo "NN"=$NN
PNTS=`echo $NN*20 | bc`
$PR 0  0  0  0 $NN $PNTS
$PR 2  0  2  0 $NN $PNTS
$PR 0  0  1  0 $NN $PNTS
$PR 0  0  2  0 $NN $PNTS
$PR 1 -1  1 -1 $NN $PNTS
$PR 1 -1  2 -1 $NN $PNTS
$PR 1  0  1  0 $NN $PNTS
$PR 1  0  2  0 $NN $PNTS
$PR 1  1  1  1 $NN $PNTS
$PR 2  1  1  1 $NN $PNTS
$PR 2 -2  2 -2 $NN $PNTS
$PR 2  1  2  1 $NN $PNTS
$PR 2  2  2  2 $NN $PNTS
rm -f RHOTMP/ro.*
rm IN_PROGRESS
touch ALL_DONE
cd ..
