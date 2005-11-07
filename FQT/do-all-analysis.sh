if [ "$1" == "" ]
then
PHIDIRS=`ls -d Phi*/`
else
PHIDIRS=`ls -d $@`
fi
if [ "$PHIDIRS"  == "" ]
then
echo "Usage: do-all-analysis.sh [Dirs_to_anaylise]"
exit 
fi
for f in $PHIDIRS
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
L=`ls Store-*-*gz`
if [ "$L" == "" ]
then
ONESTORE=`ls Store-*-*| tail -1` 
if [ "$ONESTORE" == "" ]
then
echo "[ERROR] No Store files found in" $f
exit
fi
NN=`cat $ONESTORE | awk -F : '{if ($1=="NN") print $2}'` 
else
ONESTORE=`ls Store-*-*gz| tail -1` 
#if [ "$ONESTORE" == "" ]
#then
#echo "[ERROR] No Store files found in" $f
#exit
#fi
NN=`cat $ONESTORE | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
fi
NPTS=`echo "$NN*60"| bc`
gunzip Store*gz
ls Store* | sort -t - -k 2 -k 3 -n > listamsd
$PR listamsd $NPTS
#calcola i correlatori P_n(cos(theta))
PR=$HOME/ELLIPSOIDS/FQT/calcCn
$PR listamsd $NPTS 
#gzip Store*
#gunzip Store*gz
ls Store* > junk
cat junk | while read nomefile
do
echo  $nomefile
REFTIME=`awk -F":" 'BEGIN {found=0} {if ($1=="refTime") { print $2; found=1}} END {if (!found) printf("0.0");}' $nomefile`
awk -v reft="$REFTIME" -F":" 'BEGIN {cc=0} {if (cc>0) {if ($1=="time") printf("time: %.15G\n", reft+$2); else print $0}; if ($1=="@@@") cc=cc+1; }' $nomefile > Cnf$nomefile
#awk -F":" 'BEGIN {cc=0} {if (cc>0) print $0; if ($1=="@@@") cc=cc+1; }' $nomefile > Cnf$nomefile
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
$PERC/sq-statico-1k << !
2
99
!
sh $PERC/dosqlmlpmp.sh
PR=$HOME/ELLIPSOIDS/FQT/calcfqt
echo "NN"=$NN
PNTS=`echo "$NN*20"| bc`
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
done
