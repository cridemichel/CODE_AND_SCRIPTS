PERC=$HOME/postdoc/hardellipsoid/hardellSVN/CODE/FQT/
if [ "$1" == "" ]
then
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
else
EL=$1
fi
FNT="allTAUtra-X0_${EL}.dat"
FNR="allTAUrot-X0_${EL}.dat"
echo -n "" > $FNT
echo -n "" > $FNR
A1=1.0
B1=0.1
C1=1.5
A2=1.0
B2=0.1
C2=1.5
for f in Phi* 
do
cd $f
PHI=`echo $f | awk -F Phi '{print $2}'`
if [ ! -e Fqs-10 ]
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
echo "Find Maximum of S(q)..."
MAXQ=` $PERC/findmax Sq.dat 12 | awk '{if ($1 > 29) printf("29"); else printf("%d",$1)}'`
#STA=`tail -n 50 screen_ell${EL}EQ${PHI} | awk '{if ($1=="[MSDcheck]") print $5}'` 
echo "MAXQ="$MAXQ
ST=0.00001
echo "a=$A1; b=$B1; c=$C1; fit [$ST:] a*exp(-(x/b)) \"Fqs-$MAXQ\" via a,b; fit [$ST:] a*exp(-((x/b)**c)) \"Fqs-$MAXQ\" via a,b,c; print a, b, c" > fit.tmp
echo "Fit della Fself"
gnuplot fit.tmp > gpout.tmp 2>&1  
A1=`tail -1 gpout.tmp | awk '{print $1}'`
B1=`tail -1 gpout.tmp | awk '{print $2}'`
C1=`tail -1 gpout.tmp | awk '{print $3}'`
#echo "A="$A "B="$B "C="$C
TAUM=`echo "gamma(1.0/${C1})*${B1}/${C1}" | octave | awk '{if ($1=="ans") print $3}'`
ST=0.00001
echo $EL $PHI $TAUM >> ../$FNT
echo "a=$A2; b=$B2; c=$C2; fit [$ST:] a*exp(-(x/b)) \"Cn.dat\" via a,b; fit [$ST:(3*b)] a*exp(-((x/b)**c)) \"Cn.dat\" using 1:2 via a,b,c; print a, b, c " > fit.tmp
gnuplot fit.tmp > gpout.tmp 2>&1
echo "Fit della C2"
A2=`tail -1 gpout.tmp | awk '{print $1}'`
B2=`tail -1 gpout.tmp | awk '{print $2}'`
C2=`tail -1 gpout.tmp | awk '{print $3}'`
TAUM=`echo "gamma(1.0/${C2})*${B2}/${C2}" | octave | awk '{if ($1=="ans") print $3}'`
rm gpout.tmp
echo $EL $PHI $TAUM >> ../$FNR
rm fit.tmp
cd ..
done
