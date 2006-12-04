FNT="fitall.dat"
echo -n "" > $FNT
A1=1.0
B1=1.1
C1=1.5
if [ "$1" == "" ]
then 
echo "fitallFqs.sh <lista_files>"
exit
fi
for f in `cat $1` 
do
echo "Processing file=" $f
ST=0.00001
STA=0.00001
STB=`cat $f| tail -1 | LANG=C awk '{print $1/10}'`
QVAL=`echo $f| awk -F - '{print $2+2}'`
echo "STA=" $STA "STB=" $STB "Q=" $QVAL
echo "a=$A1; b=$B1; c=$C1; fit [$STA:] a*exp(-(x/b)) \"$f\" via a,b; fit [$STA:$STB] a*exp(-((x/b)**c)) \"$f\" via a,b,c; print a, b, c" > fit.tmp
echo "Fit della Fself"
gnuplot fit.tmp > gpout.tmp 2>&1  
cp fit.tmp fitFqs-$QVAL.tmp
A1=`tail -1 gpout.tmp | awk '{print $1}'`
B1=`tail -1 gpout.tmp | awk '{print $2}'`
C1=`tail -1 gpout.tmp | awk '{print $3}'`
#echo "A="$A "B="$B "C="$C
if [ "${C1}" == "" ]
then
A1=1.0
B1=1.1
C1=1.5
echo "a=$A1; b=$B1; c=$C1; fit [$STA:] a*exp(-(x/b)) \"$f\" via a,b; fit [$STA:$STB] a*exp(-((x/b)**c)) \"$f\" via a,b,c; print a, b, c" > fit.tmp
echo "Fit della Fself"
gnuplot fit.tmp > gpout.tmp 2>&1 
cp fit.tmp fitFqs-$QVAL.tmp
A1=`tail -1 gpout.tmp | awk '{print $1}'`
B1=`tail -1 gpout.tmp | awk '{print $2}'`
C1=`tail -1 gpout.tmp | awk '{print $3}'`
fi
TAUM=`echo "gamma(1.0/${C1})*${B1}/${C1}" | octave | awk '{if ($1=="ans") print $3}'`
ST=0.00001
echo $QVAL $PHI $TAUM >> $FNT
###
rm fit.tmp
done
