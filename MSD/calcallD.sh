if [ "$1" == "" ]
then
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
else
EL=$1
fi
FNT="allDtra-X0_${EL}.dat"
FNR="allDrot-X0_${EL}.dat"
echo -n "" > $FNT
echo -n "" > $FNR
for f in Phi* 
do
cd $f
if [ ! -e rotMSDcnf.dat ]
then
echo "Phi=" $PHI " The file rotMSDcnf.dat does not exist, skipping..."
cd ..
continue
fi
if [ ! -e MSDcnf.dat ]
then
echo "Phi=" $PHI " The file MSDcnf.dat does not exist, skipping..."
cd ..
continue
fi
PHI=`echo $f | awk -F Phi '{print $2}'`
echo "Processing Phi=" $PHI
STA=`tail -n 50 screen_ell${EL}EQ${PHI} | awk '{if ($1=="[MSDcheck]") print $5}'` 
ST=`echo "$STA*10.0" | bc -l`
echo "fit [$ST:] a*x+b 'MSDcnf.dat' via a,b" > fit.tmp
gnuplot fit.tmp > gpout.tmp 2>&1  
A=`cat gpout.tmp | awk 'BEGIN {pr=0} {if (pr==1 && $1=="a") {print $3; pr++;}; if ($1=="Final") pr=1;}'`
Dt=`echo $A/6.0 | bc -l`
echo $EL $PHI $Dt >> ../$FNT
echo "fit [$ST:] a*x+b 'rotMSDcnf.dat' via a,b" > fit.tmp
gnuplot fit.tmp > gpout.tmp 2>&1
A=`cat gpout.tmp | awk 'BEGIN {pr=0} {if (pr==1 && $1=="a") {print $3; pr++;}; if ($1=="Final") pr=1;}'`
rm gpout.tmp
Dr=`echo $A/6.0 | bc -l`
echo $EL $PHI $Dr >> ../$FNR
rm fit.tmp
cd ..
done
