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
SIG="$2"
PHI=`echo $f| awk -F _ '{print $2}'`
echo "Processing file=" $f
FN="SP-popul.dat"
STA="600.0" #0.001
L=`tail -n 1 startEQ.cnf| awk '{print $1}'`
PN=`cat startEQ.cnf | grep parnum | awk -F : '{print $2}'`
echo "L=" $L " PN=" $PN
STB="3000.0" #`cat $FN| tail -1 | LANG=C awk '{print 9.*$1/10.}'`
A1=1000.0
B1=1000.1
echo "a=$A1; b=$B1; fit [${STA}:${STB}] a*exp(-(x/b)) \"$FN\" via a,b;" > fit.tmp
echo "Exp Fit"
gnuplot fit.tmp > gpout.tmp 2>&1 
A1=`tail -n 12 gpout.tmp | awk '{if ($1=="a" && $2=="=" ) print $3}'`
B1=`tail -n 12 gpout.tmp | awk '{if ($1=="b" && $2=="=" ) print $3}'`
echo "A1= " $A1 " B1= " $B1
KK=`echo "1.0/${B1}" | bc -l`
#KK=`echo "($PN/$L/$L/$L)/${B1}" | bc -l`
echo "$PHI $KK"  >> ../$FNT
###
#rm fit.tmp
cd ..
done
