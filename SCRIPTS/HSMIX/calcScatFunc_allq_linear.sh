PERC=$HOME/HSMIX/
if [ "$1" == "" ]
then
FL="lista"
else
FL="$1"
fi
if [ "$2" == "" ]
then
SKIP=""
else
SKIP="-s $2"
fi
#$PERC/calcACV $FL 1000
$PERC/calcfqtself_linear $SKIP $FL 1000 0 99 
$PERC/calcrho $FL 0 99
$PERC/calcfqtcoll_linear 0 99 1000
