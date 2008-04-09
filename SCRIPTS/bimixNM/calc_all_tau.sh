PFT=$HOME/postdoc/WOLYNES/NJP/bimixNM/find_first_below_linint
INVE=`echo "1.0/e(1.0)"|bc -l`
echo -n "" > all_tau.dat
for f in NM*-*
do
cd $f
NM=`echo $f | awk -F M '{print $2}'`
for ff in rho*
do
cd $ff
RHO=`echo $ff | awk -F o '{print $2}'`
for fff in NVE-T*
do 
cd $fff
TEMP=`echo $fff | awk -F T '{print $2}'`
#pick one, last one
FN=`ls -rt Fqs*| tail -1`
#echo "FN=" $FN
TAU=`$PFT $FN $INVE`
echo "$TEMP $TAU $RHO $NM" >> ../../../all_tau.dat
cd ..
done
cd ..
done
cd ..
done
