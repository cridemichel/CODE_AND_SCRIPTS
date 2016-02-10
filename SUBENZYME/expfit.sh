FNT="fitall.dat"
echo -n "" > $FNT
if [ "$1" == "" ]
then 
echo "expfit.sh <lista_dirs>"
exit
fi
for f in `cat $1` 
do
cd $f
PHI=`echo $f| awk -F _ '{print $2}'`
echo "Processing file=" $f
FN="SP-popul.dat"
STA=0.00001
STB=`cat $FN| tail -1 | LANG=C awk '{print 9.*$1/10.}'`
A1=1.0
B1=1.1
C1=1.5
echo "a=$A1; b=$B1; c=$C1; fit [$STA:] a*exp(-(x/b)) \"$FN\" via a,b; fit [$STA:$STB] a*exp(-((x/b)**c)) \"$FN\" via a,b,c; print a, b, c" > fit.tmp
echo "Exp Fit"
gnuplot fit.tmp > gpout.tmp 2>&1 
cp fit.tmp fitFqs-$QVAL.tmp
A1=`tail -1 gpout.tmp | awk '{print $1}'`
B1=`tail -1 gpout.tmp | awk '{print $2}'`
C1=`tail -1 gpout.tmp | awk '{print $3}'`
KK=`echo "1.0/${B1}" | bc -l`
echo "$PHI $KK" >> ../$FNT
###
rm fit.tmp
cd ..
done
