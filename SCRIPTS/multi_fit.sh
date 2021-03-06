export LANG=C
ROOTEXE="root.exe"
PERC=$HOME/postdoc/hardellipsoid/hardellSVN/CODE/SCRIPTS/
SCALFACTS="$HOME/ELLIPSOIDS/scaling_factors.dat"
GSF="$HOME/ELLIPSOIDS/get_scalfact"
GTM="$HOME/ELLIPSOIDS/getTimeFromMSD"
BEGTIME="ACV"
CNBEG=0.02
FQSBEG=0.15
FQCBEG=0.15
DOCN=1
DOFQS=0
DOFQC=0
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
if [ "$DOFQC" == "1" ]
then
echo -n "" > $FQCTAUVSX0
echo -n "" > $FQCBETAVSX0
echo -n "" > $FQCCHISQVSX0
fi
if [ "$DOFQS" == "1" ]
then
echo -n "" > $FQSTAUVSX0
echo -n "" > $FQSBETAVSX0
echo -n "" > $FQSCHISQVSX0
fi
if [ "$DOCN" == "1" ]
then
echo -n "" > $CNTAUVSX0
echo -n "" > $CNBETAVSX0
echo -n "" > $CNCHISQVSX0
fi
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
TFACT=`$GSF $SCALFACTS $X0`
if [ ! -d "Phi$2" ]
then
cd ..
continue
fi
cd Phi$2
echo "Processing " $X0 " Phi" $2 "( Scaling Factor=" $TFACT " )"
#echo $RCMD " " `pwd`
if [ \( -e Cn.dat \) -a \( "$f" != "X0_1.0" \) -a \( "$DOCN" == "1" \) ]
then
echo -n "Fitting Cn.dat..."
if [ "$BEGTIME" == "MSD" ]
then
MSDTHR="0.2"
CNBEGSC=`$GTM rotMSDcnf.dat $MSDTHR`
elif [ "$BEGTIME" == "ACV" ]
then
#CNBEGSC=`$GSF ../../tempdecad_omACV_Phi$2.dat $X0`
CNBEGSC=`$GSF ../../tempdecad_omACV_Phi$2.dat $X0`
else
CNBEGSC=`echo "$CNBEG*$TFACT"|bc -l`
fi
CNBEGSC=`echo $CNBEGSC*1.5|bc -l`
RCMD=`echo "$PERC/fitSE.C(\"Cn.dat\",$TYPE,$CNBEGSC,1000,0)"`
$ROOTEXE -b -q $RCMD > $LF
TAU=`cat $LF | tail -4 | awk -v sf=$TFACT '{if ($2=="p1") print $3/sf}'`
BETA=`cat $LF | tail -4 | awk '{if ($2=="p2") print $3}'`
CHISQ=`cat $LF | tail -4 | awk '{if ($1=="CHISQUARE:") print $2}'`
echo $X0 $TAU >> ../../$CNTAUVSX0
echo $X0 $BETA >> ../../$CNBETAVSX0
echo $X0 $CHISQ >> ../../$CNCHISQVSX0
echo "done"
fi
#
FQS=`ls -rt -1 Fqs*max 2> /dev/null | tail -1 2> /dev/null` 
if [ \( "$FQS" != "" \) -a \( "$DOFQS" == "1" \) ]
then
echo -n "Fitting " $FQS "..."
if [ "$BEGTIME" == "MSD" ]
then
MSDTHR=`echo "$TFACT/10.0" | bc -l`
FQSBEGSC=`$GTM MSDcnf.dat $MSDTHR`
elif [ "$BEGTIME" == "ACV" ]
then
FQSBEGSC=`$GSF ../../tempdecad_velACV_Phi$2.dat $X0`
else
FQSBEGSC=`echo "$FQSBEG*$TFACT"|bc -l`
fi
RCMD=`echo "$PERC/fitSE.C(\"$FQS\",$TYPE,$FQSBEGSC,1000,0)"`
$ROOTEXE -b -q $RCMD > $LF
TAU=`cat $LF | tail -4 | awk -v sf=$TFACT '{if ($2=="p1") print $3/sf}'`
BETA=`cat $LF | tail -4 | awk '{if ($2=="p2") print $3}'`
CHISQ=`cat $LF | tail -4 | awk '{if ($1=="CHISQUARE:") print $2}'`
echo $X0 $TAU >> ../../$FQSTAUVSX0
echo $X0 $BETA >> ../../$FQSBETAVSX0
echo $X0 $CHISQ >> ../../$FQSCHISQVSX0
echo "done"
fi
#
FQC=`ls -rt -1 N-sqt*max 2> /dev/null | tail -1 2> /dev/null` 
if [ "$FQC" != "" -a \( "$DOFQC" == "1" \) ]
then
echo -n "Fitting " $FQC "..."
#echo "qui1 " $FQCBEG  " " $TFACT
if [ "$BEGTIME" == "MSD" ]
then
MSDTHR=`echo "$TFACT/10.0" | bc -l`
FQCBEGSC=`$GTM MSDcnf.dat $MSDTHR`
elif [ "$BEGTIME" == "ACV" ]
then
FQCBEGSC=`$GSF ../../tempdecad_velACV_Phi$2.dat $X0`
else
FQCBEGSC=`echo "$FQCBEG*$TFACT"|bc -l`
fi
echo "FQCBEG= " $FQCBEG " FQCBEGSC= " $FQCBEGSC " SF= " $TFACT 
RCMD=`echo "$PERC/fitSE.C(\"$FQC\",$TYPE,$FQCBEGSC,1000,0)"`
$ROOTEXE -b -q $RCMD > $LF
TAU=`cat $LF | tail -4 | awk -v sf=$TFACT '{if ($2=="p1") print $3/sf}'`
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
