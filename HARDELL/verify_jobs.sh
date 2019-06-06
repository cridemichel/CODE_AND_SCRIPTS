DIR=`pwd`
AA=`basename $DIR`
X0=`echo $AA | awk -F _ '{print $2}'`
echo "X0=" $X0
NRUN=`ps ax | grep veff | grep mosrun | grep -F -e "$X0"|wc -l`
NDONE=`tail IR_*/veff*dat|grep 99900000000| wc -l`
echo "NRUN=" $NRUN
echo "NDONE=" $NDONE
echo "TOT=" $[$[NRUN]+$[NDONE]]
