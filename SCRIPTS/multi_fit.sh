ROOTEXE="root.exe"
PERC=$HOME/postdoc/hardellipsoid/hardellSVN/CODE/SCRIPTS/
CNBEG=0.1
FQSBEG=0.1
FQCBEG=0.1
DOCN=1
DOFQS=1
DOFQC=1
if [ "$3" = "" ]
then
TYPE="2"
else
TYPE=$3
fi
if [ "$TYPE" == "0" ] 
then
EXT="exp"
else
EXT="str"
fi
CNTAUVSX0="tau_vs_X0_Cn_Phi$2_$EXT.dat"
CNBETAVSX0="beta_vs_X0_Cn_Phi$2_$EXT.dat"
CNCHISQVSX0="chisq_vs_X0_Cn_Phi$2_$EXT.dat"
FQSTAUVSX0="tau_vs_X0_Fqs_Phi$2_$EXT.dat"
FQSBETAVSX0="beta_vs_X0_Fqs_Phi$2_$EXT.dat"
FQSCHISQVSX0="chisq_vs_X0_Fqs_Phi$2_$EXT.dat"
FQCTAUVSX0="tau_vs_X0_N-sqt_Phi$2_$EXT.dat"
FQCBETAVSX0="beta_vs_X0_N-sqt_Phi$2_$EXT.dat"
FQCCHISQVSX0="chisq_vs_X0_N-sqt_Phi$2_$EXT.dat"
echo -n "" > $FQCTAUVSX0
echo -n "" > $FQCBETAVSX0
echo -n "" > $FQCCHISQVSX0
echo -n "" > $FQSTAUVSX0
echo -n "" > $FQSBETAVSX0
echo -n "" > $FQSCHISQVSX0
echo -n "" > $CNTAUVSX0
echo -n "" > $CNBETAVSX0
echo -n "" > $CNCHISQVSX0
LF=_rootexe.log
if [ "$1" == "" ]
then
echo "multi_fit.sh '<list of dir>' <volume_fraction>"
exit
fi
for f in $1
do
cd $f
X0=`echo $f | awk -F '_' '{print $2}'`
if [ ! -d "Phi$2" ]
then
cd ..
continue
fi
cd Phi$2
echo "Processing " $X0 " Phi" $2 
#echo $RCMD " " `pwd`
if [ \( -e Cn.dat \) -a \( "$f" != "X0_1.0" \) -a \( "$DOCN" == "1" \) ]
then
echo -n "Fitting Cn.dat..."
RCMD=`echo "$PERC/fitSE.C(\"Cn.dat\",$TYPE,$CNBEG,1000,0)"`
$ROOTEXE -b -q $RCMD > $LF
TAU=`cat $LF | tail -4 | awk '{if ($2=="p1") print $3}'`
BETA=`cat $LF | tail -4 | awk '{if ($2=="p2") print $3}'`
CHISQ=`cat $LF | tail -4 | awk '{if ($1=="CHISQUARE:") print $2}'`
echo $X0 $TAU >> ../../$CNTAUVSX0
echo $X0 $BETA >> ../../$CNBETAVSX0
echo $X0 $CHISQ >> ../../$CNCHISQVSX0
echo "done"
fi
#
FQS=`ls -r -1 Fqs*max 2> /dev/null | tail -1 2> /dev/null` 
if [ \( "$FQS" != "" \) -a \( "$DOFQS" == "1" \) ]
then
echo -n "Fitting " $FQS "..."
RCMD=`echo "$PERC/fitSE.C(\"$FQS\",$TYPE,$FQSBEG,1000,0)"`
$ROOTEXE -b -q $RCMD > $LF
TAU=`cat $LF | tail -4 | awk '{if ($2=="p1") print $3}'`
BETA=`cat $LF | tail -4 | awk '{if ($2=="p2") print $3}'`
CHISQ=`cat $LF | tail -4 | awk '{if ($1=="CHISQUARE:") print $2}'`
echo $X0 $TAU >> ../../$FQSTAUVSX0
echo $X0 $BETA >> ../../$FQSBETAVSX0
echo $X0 $CHISQ >> ../../$FQSCHISQVSX0
echo "done"
fi
#
FQC=`ls -r -1 N-sqt*max 2> /dev/null | tail -1 2> /dev/null` 
if [ "$FQC" != "" -a \( "$DOFQC" == "1" \) ]
then
echo -n "Fitting " $FQC "..."
RCMD=`echo "$PERC/fitSE.C(\"$FQC\",$TYPE,$FQCBEG,1000,0)"`
$ROOTEXE -b -q $RCMD > $LF
TAU=`cat $LF | tail -4 | awk '{if ($2=="p1") print $3}'`
BETA=`cat $LF | tail -4 | awk '{if ($2=="p2") print $3}'`
CHISQ=`cat $LF | tail -4 | awk '{if ($1=="CHISQUARE:") print $2}'`
echo $X0 $TAU >> ../../$FQCTAUVSX0
echo $X0 $BETA >> ../../$FQCBETAVSX0
echo $X0 $CHISQ >> ../../$FQCCHISQVSX0
echo "done"
fi
if [ -e $LF ] 
then
rm $LF
fi
cd ..
cd ..
done
