#!/bin/bash
#normalmente se i file ci sono già non ricalcola le varie grandezze ma si puo' forzare il ricalcolo
FORCE_MSD=0
FORCE_FQSELF=0
FORCE_SQ=0
FORCE_CN=0
EXE_PATH=$HOME/ELLIPSOIDS/FQT/
QMIN=2
QMAX=30
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
if [ -e IGNORE_THIS ]
then 
cd ..
continue
fi
if [ -e ALL_DONE_SELF ]
then 
cd ..
continue
fi
touch IN_PROGRESS
PR=$EXE_PATH/MSD/calcmsd
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
#if [ "$ONESTORE" == "" ]
#then
#echo "[ERROR] No Store files found in" $f
#exit
#fi
NN=`cat $ONESTORE | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
fi
NPTS=`echo "$NN*60"| bc`
if [ "$L" != "" ]
then
gunzip Store*gz
fi
ls Store* | sort -t - -k 2 -k 3 -n > listamsd
if [ ! \( -e MSDcnf.dat \) -o \( $FORCE_MSD == "1" \) ]
then
echo "CALCOLO MSD..."
$PR listamsd $NPTS
echo "DONE"
fi
if [ ! \( -e Cn.dat  \) -o \( $FORCE_CN == "1" \) ]
then
#calcola i correlatori P_n(cos(theta))
PR=$EXE_PATH/FQT/calcCn
echo "CALCOLO Cn..."
$PR listamsd $NPTS 
echo "DONE"
fi
if [ ! \( -e Sq.dat \) -o \( $FORCE_SQ == "1" \) ]
then
echo "CALCOLO SQ..."
ls Store-*-0 > listasq
PR=$EXE_PATH/FQT/calcSq
$PR --cnf listasq
echo "DONE"
fi
if [ ! \( -e Fqs-10 \) -o \( $FORCE_FQSELF == "1" \) ]
then
PR=$EXE_PATH/FQT/calcfqtself
echo "CALCOLO Fself..."
$PR listamsd $NPTS $QMIN $QMAX
echo "DONE"
fi
# =====================================================
if [ "$L" != "" ] 
then
gzip Store-*-*
fi
rm IN_PROGRESS
touch ALL_DONE
cd ..
done
