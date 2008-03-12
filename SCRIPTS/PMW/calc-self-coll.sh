#!/bin/bash
#normalmente se i file ci sono gi� non ricalcola le varie grandezze ma si puo' forzare il ricalcolo
FORCE_MSD=1
FORCE_FQSELF=0
FORCE_FQCOLL=0
FORCE_SQ=0
FORCE=0
if [ $FORCE_FQSELF == "1" ]
then 
FORCE=1
fi
if [ $FORCE_FQCOLL == "1" ]
then 
FORCE=1
fi
if [ $FORCE_MSD == "1" ]
then 
FORCE=1
fi
if [ $FORCE_SQ == "1" ]
then 
FORCE=1
fi
EXE_PATH=$HOME/MD/PWM/
QMIN=0
QMAX="-1"
if [ "$1" == "" ]
then
PHIDIRS=`ls -d Phi*/`
else
PHIDIRS=`ls -d $@`
fi
if [ "$PHIDIRS"  == "" ]
then
echo "Usage: calc-self.sh [Dirs_to_anaylise]"
exit 
fi
for f in $PHIDIRS
do
cd $f
if [ $FORCE == "0" ]
then
if [ \( -e Sq.dat \) -a \( -e Fqs-0 \) -a \( -e MSDcnf.dat \) -a \( N-sqt.k=000 \)]
then
echo "all done, skipping Phi=" $f
cd ..
continue
fi
fi
if [ -e IGNORE_THIS ]
then 
cd ..
continue
fi
#if [ -e ALL_DONE_SELF ]
#then 
#cd ..
#continue
#fi
echo "Doing Phi=" $f
touch IN_PROGRESS
PR=$EXE_PATH/calcmsd
L=`ls Store-*-0.gz 2> /dev/null`
#echo "L=" $L
if [ "$L" == "" ]
then
ONESTORE=`ls Store-*-0| tail -1` 
if [ "$ONESTORE" == "" ]
then
echo "[ERROR] No Store files found in" $f
exit
fi
NN=`cat $ONESTORE | awk -F : '{if ($1=="NN") print $2}'` 
else
ONESTORE=`ls Store-*-0.gz| tail -1` 
#if [ "$ONESTORE" == "" ]
#then
#echo "[ERROR] No Store files found in" $f
#exit
#fi
NN=`cat $ONESTORE | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
fi
NPTS=`echo "$NN*60"| bc`
#echo "NPTS=" $NPTS
#exit
#if [ "$L" != "" ]
#then
#gunzip -f Store*gz
#fi
rm -f listamsd
i=0
while [ 1 -eq 1 ]
do
if [ ! -e Store-$i-0.gz -a ! -e Store-$i-0 ]
then 
break
fi
if [ "$L" != "" ]
then
gunzip -f Store-$i-*gz
fi
ls Store-$i-* | sort -t - -k 2 -k 3 -n >> listamsd
i=$[$i+1]
done
if [ ! \( -e MSDcnf.dat \) -o \( $FORCE_MSD == "1" \) ]
then
echo "CALCOLO MSD..."
$PR listamsd $NPTS
echo "DONE"
fi
if [ ! \( -e Sq.dat \) -o \( $FORCE_SQ == "1" \) ]
then
echo "CALCOLO SQ..."
ls Store-*-0 > listasq
PR=$EXE_PATH/calcSq
$PR --cnf listasq
echo "DONE"
fi
if [ "$QMAX" == "-1" ]
then
PR=$EXE_PATH/findmax
SQMAX=`$PR Sq.dat 0`
SQMAXFQS=`echo "$SQMAX" | awk -v sqmax=$SQMAX '{printf("%d",sqmax)}'`
echo "QMAX FROM SQ= "$SQMAX " [" $SQMAXFQS "]"
fi
if [ ! \( -e Fqs-0 \) -o \( $FORCE_FQSELF == "1" \) ]
then
PR=$EXE_PATH/calcfqtself
echo "CALCOLO Fself..."
if [ "$QMAX" == "-1" ]
then
$PR listamsd $NPTS 0 0
$PR listamsd $NPTS $SQMAX $SQMAX
cp Fqs-$SQMAXFQS Fqs-$SQMAXFQS.max
else
$PR listamsd $NPTS $QMIN $QMAX
fi
echo "DONE"
fi
if [ ! \( -e N-sqt.k=000 \) -o \( $FORCE_FQCOLL == "1" \) ]
then
PR=$EXE_PATH/calcrho
PR2=$EXE_PATH/calcfqtcoll
echo "CALCOLO Fcoll..."
if [ "$QMAX" == "-1" ]
then
$PR listamsd 0 0 
$PR listamsd $SQMAX $SQMAX
$PR2 --ncomps 1 0 0 $NPTS
$PR2 --ncomps 1 $SQMAX $SQMAX $NPTS
rm -fr RHOTMP/
cp N-sqt.k=$SQMAX N-sqt.k=$SQMAX.max
else
$PR listamsd $QMIN $QMAX
$PR2 --ncomps 1 $QMIN $QMAX $NPTS
rm -fr RHOTMP/
fi
echo "DONE"
fi
# =====================================================
i=0
while [ 1 -eq 1 ]
do
if [ ! -e Store-$i-0 ]
then 
break
fi
if [ "$L" != "" ]
then
gzip -f Store-$i-*
fi
i=$[$i+1]
done
#if [ "$L" != "" ] 
#then
#gzip -f Store-*-*
#fi
rm IN_PROGRESS
#touch ALL_DONE
cd ..
done
