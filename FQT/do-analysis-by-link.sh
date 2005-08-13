SCR=../FQT/do-all-analysis.sh
if [ "$1" == "" ]
then
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
else
EL=$1
fi
SCRLNK=do-analysis-X0_${EL}.sh
ln -sf $SCR $SCRLNK
sh $SCRLNK > screen_$SCRLNK 2>&1 &
