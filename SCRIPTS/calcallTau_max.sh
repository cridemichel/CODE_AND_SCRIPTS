PERC=$HOME/ELLIPSOIDS/FQT/
if [ "$2" == "" ]
then
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
else
EL=$2
fi
FNTM="allTAUtraM-X0_${EL}.dat"
FNT="allTAUtra-X0_${EL}.dat"
FNR="allTAUrot-X0_${EL}.dat"
FNTMS="allTAUtraM-X0_${EL}_2D.dat"
FNTS="allTAUtra-X0_${EL}_2D.dat"
FNRS="allTAUrot-X0_${EL}_2D.dat"
echo -n "" > $FNT
echo -n "" > $FNTM
echo -n "" > $FNR
echo -n "" > $FNTS
echo -n "" > $FNTMS
echo -n "" > $FNRS
A1=1.0
B1=1.1
C1=1.5
A3=1.0
B3=1.1
C3=1.5
A2=1.0
B2=1.1
C2=1.5
if [ "$1" == "" ]
then
DIRS='Phi*'
else
DIRS="$1"
fi
#echo "DIRS:" $DIRS
for f in $DIRS
do
cd $f
PHI=`echo $f | awk -F Phi '{print $2}'`
if [ ! -e Fqs-0 ]
then
echo "Phi=" $PHI " The file Fqs-10 does not exist, skipping..."
cd ..
continue
fi
if [ ! -e Cn.dat ]
then
echo "Phi=" $PHI " The file Cn.dat does not exist, skipping..."
cd ..
continue
fi
echo "Processing Phi=" $PHI
#echo "Find Maximum of S(q)..."
#MAXQF=` $PERC/findmax Sq.dat 2 | awk '{if ($1 > 29) printf("29"); else printf("%d",$1)}'`
#MAXQA=$[$MAXQF-1]
#MAXQB=$[$MAXQF+1]
#GPFILE="fitsq.tmp"
#echo "maxqa=$MAXQA; maxqb=$MAXQB; if (maxqa < 2) maxqa=2; else; a=1; b=1; c=1; fit [maxqa:maxqb] a*x*x+b*x+c \"Sq.dat\" via a,b,c; maxq=-b/(2*a); if (maxq-floor(maxq)>=0.5) print ceil(maxq); else print floor(maxq)" > $GPFILE
#gnuplot $GPFILE > gpoutsq.tmp 2>&1
#MAXQ=`tail -1 gpoutsq.tmp| awk -v maxqa=$MAXQA -v maxqb=$MAXQB -v maxqf=$MAXQF '{if ($1=="" || $1 > maxqb || $1 < maxqa) printf("%d", maxqf); else printf("%d",$1)}'`
#MAXQ2=$MAXQ
#MAXQ=`echo $MAXQ2 | awk -v maxq=$MAXQ2 '{if (maxq > 29) print "29"; else if (maxq < 2) print "2"; else print maxq}'`
#STA=`tail -n 50 screen_ell${EL}EQ${PHI} | awk '{if ($1=="[MSDcheck]") print $5}'` 
MAXQ=`ls Fqs-*max | awk -F . '{print $1}'| awk -F - '{print $2}'`
echo "MAXQ="$MAXQ
ST=0.00001
STA=0.00001
STB=`cat Fqs-$MAXQ.max| tail -1 | LANG=C awk '{print $1/10}'`
#MAXQ=2
echo "STA=" $STA "STB=" $STB
echo "a=$A1; b=$B1; c=$C1; fit [$STA:] a*exp(-(x/b)) \"Fqs-$MAXQ.max\" via a,b; fit [$STA:$STB] a*exp(-((x/b)**c)) \"Fqs-$MAXQ.max\" via a,b,c; print a, b, c" > fit.tmp
echo "Fit della Fself"
gnuplot fit.tmp > gpout.tmp 2>&1  
cp fit.tmp fitFqs-$MAXQ.tmp
A1=`tail -1 gpout.tmp | awk '{print $1}'`
B1=`tail -1 gpout.tmp | awk '{print $2}'`
C1=`tail -1 gpout.tmp | awk '{print $3}'`
#echo "A="$A "B="$B "C="$C
if [ "${C1}" == "" ]
then
A1=1.0
B1=1.1
C1=1.5
echo "a=$A1; b=$B1; c=$C1; fit [$STA:] a*exp(-(x/b)) \"Fqs-$MAXQ.max\" via a,b; fit [$STA:$STB] a*exp(-((x/b)**c)) \"Fqs-$MAXQ.max\" via a,b,c; print a, b, c" > fit.tmp
echo "Fit della Fself"
gnuplot fit.tmp > gpout.tmp 2>&1 
cp fit.tmp fitFqs-$MAXQ.tmp
A1=`tail -1 gpout.tmp | awk '{print $1}'`
B1=`tail -1 gpout.tmp | awk '{print $2}'`
C1=`tail -1 gpout.tmp | awk '{print $3}'`
fi
TAUM=`echo "gamma(1.0/${C1})*${B1}/${C1}" | octave | awk '{if ($1=="ans") print $3}'`
ST=0.00001
echo $EL $PHI $TAUM >> ../$FNTM
echo $PHI $TAUM >> ../$FNTMS
###
MAXQ=0
echo "MAXQ="$MAXQ
ST=0.00001
STA=0.00001
STB=`cat Fqs-$MAXQ| tail -1 | LANG=C awk '{print $1/10}'`
#MAXQ=2
echo "STA=" $STA "STB" $STB
echo "a=$A3; b=$B3; c=$C3; fit [$STA:] a*exp(-(x/b)) \"Fqs-$MAXQ\" via a,b; fit [$STA:$STB] a*exp(-((x/b)**c)) \"Fqs-$MAXQ\" via a,b,c; print a, b, c" > fit.tmp
echo "Fit della Fself"
gnuplot fit.tmp > gpout.tmp 2>&1  
cp fit.tmp fitFqs-0.tmp
A3=`tail -1 gpout.tmp | awk '{print $1}'`
B3=`tail -1 gpout.tmp | awk '{print $2}'`
C3=`tail -1 gpout.tmp | awk '{print $3}'`
if [ "${C3}" == "" ]
then
A3=1.0
B3=1.1
C3=1.5
echo "a=$A3; b=$B3; c=$C3; fit [$STA:] a*exp(-(x/b)) \"Fqs-$MAXQ\" via a,b; fit [$STA:$STB] a*exp(-((x/b)**c)) \"Fqs-$MAXQ\" via a,b,c; print a, b, c" > fit.tmp
echo "Fit della Fself"
gnuplot fit.tmp > gpout.tmp 2>&1  
cp fit.tmp fitFqs-0.tmp
A3=`tail -1 gpout.tmp | awk '{print $1}'`
B3=`tail -1 gpout.tmp | awk '{print $2}'`
C3=`tail -1 gpout.tmp | awk '{print $3}'`
fi
#echo "A="$A "B="$B "C="$C
TAUM=`echo "gamma(1.0/${C3})*${B3}/${C3}" | octave | awk '{if ($1=="ans") print $3}'`
ST=0.00001
echo $EL $PHI $TAUM >> ../$FNT
echo $PHI $TAUM >> ../$FNTS
###
echo "a=$A2; b=$B2; c=$C2; fit [$ST:] a*exp(-(x/b)) \"Cn.dat\" via a,b; fit [$ST:(3.0*b)] a*exp(-((x/b)**c)) \"Cn.dat\" using 1:2 via a,b,c; print a, b, c " > fit.tmp
gnuplot fit.tmp > gpout.tmp 2>&1
cp fit.tmp fitC2.tmp
echo "Fit della C2"
A2=`tail -1 gpout.tmp | awk '{print $1}'`
B2=`tail -1 gpout.tmp | awk '{print $2}'`
C2=`tail -1 gpout.tmp | awk '{print $3}'`
if [ "${C2}" == "" ]
then
A2=1.0
B2=1.1
C2=1.5
echo "a=$A2; b=$B2; c=$C2; fit [$ST:] a*exp(-(x/b)) \"Cn.dat\" via a,b; fit [$ST:(3.0*b)] a*exp(-((x/b)**c)) \"Cn.dat\" using 1:2 via a,b,c; print a, b, c " > fit.tmp
gnuplot fit.tmp > gpout.tmp 2>&1
cp fit.tmp fitC2.tmp
echo "Fit della C2"
A2=`tail -1 gpout.tmp | awk '{print $1}'`
B2=`tail -1 gpout.tmp | awk '{print $2}'`
C2=`tail -1 gpout.tmp | awk '{print $3}'`
fi
TAUM=`echo "gamma(1.0/${C2})*${B2}/${C2}" | octave | awk '{if ($1=="ans") print $3}'`
##rm gpout.tmp
echo $EL $PHI $TAUM >> ../$FNR
echo $PHI $TAUM >> ../$FNRS
rm fit.tmp
cd ..
done
