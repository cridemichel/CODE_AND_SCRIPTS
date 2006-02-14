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
MAXQ=` $PERC/findmax Sq.dat 12 | awk '{printf("%d",$1)}'`
#STA=`tail -n 50 screen_ell${EL}EQ${PHI} | awk '{if ($1=="[MSDcheck]") print $5}'` 
echo "MAXQ="$MAXQ
ST=0.1
echo "a=1.0; b=0.1; c=1.0; fit [$ST:] a*exp(-((x/b)**c)) 'Fqs-$MAXQ' via a,b,c" > fit.tmp
echo "Fit della Fself"
gnuplot fit.tmp > gpout.tmp 2>&1  
B=`cat gpout.tmp | awk 'BEGIN {pr=0} {if (pr==1 && $1=="b") {print $3; pr++;}; if ($1=="Final") pr=1;}'`
C=`cat gpout.tmp | awk 'BEGIN {pr=0} {if (pr==1 && $1=="c") {print $3; pr++;}; if ($1=="Final") pr=1;}'`
#echo "A="$A "B="$B "C="$C
TAUM=`echo "gamma(${C})*${B}/${C}" | octave | awk '{if ($1=="ans") print $3}'`
echo $EL $PHI $TAUM >> ../$FNT
echo "a=1.0; b=0.1; c=1.0; fit [$ST:] a*exp(-((x/b)**c)) 'Cn.dat' using 1:2 via a,b,c" > fit.tmp
gnuplot fit.tmp > gpout.tmp 2>&1
echo "Fit della C2"
B=`cat gpout.tmp | awk 'BEGIN {pr=0} {if (pr==1 && $1=="b") {print $3; pr++;}; if ($1=="Final") pr=1;}'`
C=`cat gpout.tmp | awk 'BEGIN {pr=0} {if (pr==1 && $1=="c") {print $3; pr++;}; if ($1=="Final") pr=1;}'`
TAUM=`echo "gamma(${C})*${B}/${C}" | octave | awk '{if ($1=="ans") print $3}'`
rm gpout.tmp
echo $EL $PHI $TAUM >> ../$FNR
rm fit.tmp
cd ..
done
