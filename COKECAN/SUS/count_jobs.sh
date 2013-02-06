NJ=`ps ax | grep mosrun | grep HC_N| wc -l`
#ND=`ls -d N_*`
echo "#jobs =" $NJ
if [ "$1" == "1" ]
then
ND=`ls -d N_*| wc -l`
CHECK=`echo "" | awk -v nd=$ND -v nj=$NJ '{if (nj < nd) print "1"; else print "0";}'`
echo "CHECK=" $CHECK
if [ "$CHECK" == "1" ]
then
echo "List of simulations not running:"
ps ax | grep mosrun > _aaa_
for f in `ls -d N_*`
do
EXIST=`cat _aaa_ |grep $f | awk '{print $1}'`
rm _aaa_
#echo "EXIST= " $EXIST
if [ "$EXIST" == "" ]
then
echo "run in dir " $f " is missing"
fi
done
fi
fi
