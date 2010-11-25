PERC=$HOME/HSMIX/
if [ "$1" == "" ]
then
FL="lista"
else
FL="$1"
fi
if [ "$2"=="" ]
then
SKIP=""
else
SKIP="-s $2"
fi
$PERC/calcfqtself $SKIP $FL 1000 0 99 
$PERC/calcrho $FL 0 99
$PERC/calcfqtcoll 0 99 1000
