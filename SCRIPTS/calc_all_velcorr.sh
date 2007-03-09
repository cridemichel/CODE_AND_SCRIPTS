if [ "$1" == "" ]
then
echo "calc_all_velcorr.sh <volume_fraction>"
fi
for f in X0_*
do 
cd $f 
if [ ! -d Phi$1 ]
then
cd ..
continue
fi
cd Phi$1
echo "Doing " $f
touch IN_PROGRESS
PR=$EXE_PATH/AUTOCORR_VEL/calcACV
L=`ls Store-*-*gz 2> /dev/null`
#echo "L=" $L
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
NN=`cat $ONESTORE | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
fi
NPTS=`echo "$NN*10"| bc`
#echo "NPTS=" $NPTS
#exit
if [ "$L" != "" ]
then
gunzip -f Store*gz
fi
ls Store-*-* | sort -t - -k 2 -k 3 -n > lista
$PR lista $NPTS
rm lista
cd ..
cd ..
done
