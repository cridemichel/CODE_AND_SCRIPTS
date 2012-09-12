# $1 =pressure
# $2 = cicli di produzione (1 ciclo = 1 tempo di equilibratura) (se < 0 non fa la produzione)
# $3 = elongazione
# $4 = temperatura
FIXVB="1"
if [ "$1" = "" ]
then
echo "Syntax: sim1statepnt <pressure> <steps> [elongation] [temperature]"
exit 0
fi
if [ "$2" = "" ]
then
echo "You have to supply the number of steps"
exit 0
fi
if [ "$3" = "" ]
then 
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
#echo "You have to supply an elongation!"
#echo "Elongazione: $EL"
#exit 0
else
EL=$3
fi
if [ "$4" == "" ]
then
TEMP="0.15"
else
TEMP="$4"
fi
if [ "$5" == "" ]
then
SIGMA="0.9"
else
SIGMA="$5"
fi
if [ "$6" != "" ]
then
FIXVB="$6"
fi
if [ "$FIXVB" == "0" ]
then
DIRSIM="sigma_${SIGMA}_noFixVb_Press$1_T$TEMP"
elif [ "$FIXVB" == "2" ]
then
DIRSIM="sigma_${SIGMA}_modB_Press$1_T$TEMP"
if [ "$SIGMA" == "-1" ]
then
SIGMA="0.608333"
fi
else
DIRSIM="sigma_${SIGMA}_Press$1_T$TEMP"
fi
PARFILE="ellipsoid_flex.par"
if [ ! -e $DIRSIM ]
then 
mkdir $DIRSIM
fi
cp $PARFILE $DIRSIM
cd $DIRSIM
rm -f COORD_TMP*
rm -f Store-*
PRESS="$1"
ELLEXE="../../ellipsoid"
INITEMP="2.0"
SIMPR="ellMC-NPT${EL}PR$1T${TEMP}SIG$SIGMA"
MOSRUN="mosrun"
#per ora il salvataggio è lineare
#=========== >>> PARAMETRI <<< =============
STORERATE="50.0"
USENNL=1
INIFILE="start.cnf"
PARNUM=1000
DT="0.05"
NNPR=1
BASEPR=1
INTSUM="20.0"
OCTAVE="octave"
OCTFILE="script.oct"
#============================================
#N.B. it's supposed that we use NNL here!!
PROL=`echo $EL | awk '{if ($0 >= 1.0) printf("1"); else printf("0");}'`
if [ $PROL -eq 1 ]
then
A0=$EL
B0=1.0
C0=1.0
RNNL=0.3
#INIL=`echo "2.0*e(1.0/3.0*l($PARNUM))*$A0" | bc -l`
#echo "qui INIL=" $INIL
if [ $USENNL -eq 0 ]
then
RCUT=`echo "2.0*$A0*1.01" | bc -l` 
fi
else
A0="1.0"
B0="$EL" 
C0="$EL"
RNNL=0.3
#INIL=`echo "5.0*e(1.0/3.0*l($PARNUM))*$B0" | bc -l`
if [ $USENNL -eq 0 ]
then
RCUT=`echo "2.0*$B0*1.01" | bc -l` 
NNLPAR="0"
fi
fi
cp ../$INIFILE .
#============ >>> set particles number and elongation in $INIFILE <<< ==================
cat $INIFILE | awk -v np=$PARNUM -v a=$A0 -v b=$B0  -v c=$C0 ' BEGIN {NAT=0} {if (NAT==1 && NR==NL+1) {print np;} else if (NAT==1 && NR==NL+2) { print (a,b,c);} else {print $0;} if ($0=="@@@") {NAT+=1; NL=NR} }' > _aaa_
mv _aaa_ $INIFILE
#============ >>> adjust spot according to $SIGMA <<< =============
if [ "$6" != "" ]
then
VB="0.00706858347057703"
echo 'function y=sphCap(h)' > $OCTFILE
echo 'y=2*acos(0)*h^2*(3*'$SIGMA'/2-h)/3-'$VB';' >> $OCTFILE
echo "endfunction" >> $OCTFILE
echo '[h, info]'"= fsolve (\"sphCap\",""$SIGMA"'/4)' >> $OCTFILE
echo "printf(\"hh %.8G\n\",h);" >> $OCTFILE
$OCTAVE -q $OCTFILE > _aaa_
SIG2=`echo $SIGMA | awk '{print $1/2}'`
HVAL=`cat _aaa_ | awk '{if ($1=="hh") print $2}'` 
ISNEG=`echo $HVAL |awk '{if ($1 <= 0) print "1"; else print "0";}'`
ISBIG=`echo $HVAL |awk -v sig2=$SIG2 '{if ($1 > sig2) print "1"; else print "2";}'`
if [ "$ISNEG" == "1" ]
then
echo "HVAL is negative ( HVAL=" $HVAL " ) exiting!"
exit
fi
if [ "$ISBIG" == "1" ]
then
echo "HVAL is too big ( HVAL=" $HVAL " ) exiting!"
exit
fi
if [ "$FIXVB" == "0" ]
then
DEL="$EL"
elif [ "$FIXVB" == "2" ]
then
DEL=`echo "0.15"|awk -v el=$EL -v sig=$SIGMA '{printf("%.8G",el-(sig/2-$1));}'`
else
DEL=`echo $HVAL | awk -v el=$EL -v sig=$SIGMA '{printf("%.8G",el-(sig/2-$1));}'`
fi
echo "DEL= " $DEL " HVAL= " $HVAL
cat $INIFILE | awk -v del=$DEL -v sig=$SIGMA -v el=$EL 'BEGIN {NAT=0} {if (NAT==1 && NR==NL+6) {printf("%.8G 0.0 0.0 %.8G\n", del, sig);} else if (NAT==1 && NR==NL+7) { printf("%.8G 0.0 0.0 %.8G\n", -del, sig);} else {print $0;}; if ($0=="@@@") {NAT+=1; NL=NR}}' > _aaa_
mv _aaa_ $INIFILE
fi
#==================================================================
if [ $USENNL -eq 1 ]
then
#usa le NNL con sticky spots!
NNLPAR="3"
#RCUT=`echo "2.0*e(0.5*l(($A0+$RNNL)*($A0+$RNNL)+($B0+$RNNL)*($B0+$RNNL)+($C0+$RNNL)*($C0+$RNNL)))" | bc -l`
RCUT="-1"
fi
echo "RCUT=" $RCUT " " "A=" $A0 "B=" $B0 "C=" $C0 "RNNL=" $RNNL "EL=" $EL
#RANDOMIZZAZIONE INIZIALE
#cp $PARFILE rand_$PARFILE
rm -f $OCTFILE
../set_params.sh $PARFILE temperat $TEMP P $PRESS bakStepAscii 5000 stepnum $2 inifile $INIFILE endfile ${SIMPR}.cor
ln -sf $ELLEXE $SIMPR
$MOSRUN ./$SIMPR -fa ./$PARFILE > screen_$SIMPR 
cd ..
