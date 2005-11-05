# $1 = volume fraction
# $2 = cicli di produzione (1 ciclo = 1 tempo di equilibratura) (se < 0 non fa la produzione)
# $3 = elongazione
# $4 = se < 0 non fa l'equilibratura, se = 0 usa il MSD per terminare l'equilibratura, se > 0 fa 
#      $4 passi di equilibratura ( 0 = default ) 
if [ "$1" = "" ]
then
echo "Syntax: production_run_extend.sh <Volume_Fraction> <factor> [elongation]"
exit 0
fi
if [ "$2" = "" ]
then
echo "You have to supply the number of production cycles"
exit 0
fi
if [ "$4" = "" ]
then 
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
#echo "You have to supply an elongation!"
#echo "Elongazione: $EL"
#exit 0
else
EL=$4
fi
cd Phi$1-PR
STOREFILE=`ls Store-*-* | sort -n -t - -k 2 -k 3 | tail -1`
GZIPPED=`echo $STOREFILE | awk -F . '{print $2}'`
if [ "$GZIPPED" == "" ]
then
cp $STOREFILE ellips.cnf
else
cp $STOREFILE ellips.cnf.gz
gunzip ellips.cnf.gz
fi
rm -f COORD_TMP*
ELLEXE="../../ellipsoid"
SIMPR="ell${EL}PR$1"
STORERATE="0.01"
NN=`echo `
#RANDOMIZZAZIONE INIZIALE
#>>> SET TEMPERATURE TO 1.0
STPS=`cat ellips.cnf | awk -F : '{ if ($1=="totStep") print $2}'`
XSTPS=`echo "$STPS $2" | awk '{printf("%d",$1*$2)}'`
cat ellips.cnf | awk -v ts=$XSTPS -F : '{if ($1=="totStep") printf("totStep: %d\n", ts); else print $0}' > junk
mv junk ellips.cnf
$SIMPR -ca ellips.cnf >> screen_$SIMPR 
cd ..
