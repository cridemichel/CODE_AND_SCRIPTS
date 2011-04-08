if [ "$1" == "" ]
then
echo "calc_cov_s.sh <s'>"
exit
fi
fn="cov_vs_s_$1.dat"
echo -n "" > $fn
for f in $1-?
do
cd $f
COV=`cat covolume.dat | tail -1 | awk '{print $2}'`
S=`echo $f| awk -F - '{print $2}'`
echo $S $COV >> ../$fn
cd ..
done
