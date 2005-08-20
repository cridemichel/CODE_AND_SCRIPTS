SCR=../FQT/do-one-fqt.sh
if [ "$1" == "" ]
then
echo "You must supply at lease the volume fraction"
exit
fi
if [ "$2" == "" ]
then
CDIR=`pwd`
DIR=`basename $CDIR`
EL=`echo $DIR | awk -F _ '{print $2}'`
else
EL=$2
fi
PHI=$1
SCRLNK=do-sqt-X0_${EL}Phi$1.sh
ln -sf $SCR $SCRLNK
sh $SCRLNK $1 > screen_$SCRLNK 2>&1 &
