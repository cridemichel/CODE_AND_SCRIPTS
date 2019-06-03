NRUN=`ps ax | grep veff | grep mosrun | grep -e '1\.3'|wc -l`
NDONE=`tail IR_*/veff*dat|grep 99900000000| wc -l`
echo "NRUN=" $NRUN
echo "NDONE=" $NDONE
echo "TOT=" $[$[NRUN]+$[NDONE]]
