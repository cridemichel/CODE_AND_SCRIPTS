DIR=`pwd`
AA=`basename $DIR`
X0=`echo $AA | awk -F _ '{print $2}'`
echo "X0=" $X0
ps ax > _aaa_
NRUN=`cat _aaa_ | grep veff | grep mosrun | grep -F -e "$X0"|wc -l`
rm _aaa_
NDONE=`tail IR_*/veff*dat|grep 99900000000| wc -l`
echo "NRUN=" $NRUN
echo "NDONE=" $NDONE
echo "TOT=" $[$[NRUN]+$[NDONE]]
