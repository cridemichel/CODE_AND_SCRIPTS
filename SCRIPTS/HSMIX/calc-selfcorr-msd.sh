for f in Phi*
do
cd $f
if [ -e IGNORE_THIS ]
then 
cd ..
continue
fi
if [ -e ALL_DONE ]
then 
cd ..
continue
fi
touch IN_PROGRESS
PR=$HOME/ELLIPSOIDS/MSD/calcmsd
if [ "$1" == "" ]
then
NN=`cat Store-0-0.gz | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
else
NN=$1
fi
NN=`echo "$NN*30"| bc`
gunzip Store*gz
ls Store* | sort -t - -k 2 -k 3 -n > listamsd
$PR listamsd $NN 
PR=$HOME/ELLIPSOIDS/FQT/calcCn
$PR listamsd $NN
#gzip Store*
#gunzip Store*gz
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
$PERC/rho-dip-0000 << !
$q
!
q=$[$q+1]
done
rm -f CnfStore* 
$PERC/sq-statico-1k << !
2
99
!
sh $PERC/dosqlmlpmp.sh
PR=$HOME/ELLIPSOIDS/FQT/calcfqt
if [ "$1" == "" ]
then
NN=`cat Store-0-0.gz | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
else
NN=$1
fi
echo "NN"=$NN
PNTS=`echo "$NN*20"| bc`
$PR 0  0  0  0 $NN $PNTS 
#$PR 2  0  2  0 $NN $PNTS
#$PR 0  0  1  0 $NN $PNTS
#$PR 0  0  2  0 $NN $PNTS
#$PR 1 -1  1 -1 $NN $PNTS
#$PR 1 -1  2 -1 $NN $PNTS
#$PR 1  0  1  0 $NN $PNTS
#$PR 1  0  2  0 $NN $PNTS
#$PR 1  1  1  1 $NN $PNTS
#$PR 2  1  1  1 $NN $PNTS
#$PR 2 -2  2 -2 $NN $PNTS
#$PR 2  1  2  1 $NN $PNTS
#$PR 2  2  2  2 $NN $PNTS
rm -f RHOTMP/ro.*
rm IN_PROGRESS
touch ALL_DONE
cd ..
done
