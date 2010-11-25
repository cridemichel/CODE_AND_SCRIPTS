MAX=$1
PERC=$HOME/HSMIX/
if [ "$2" == "" ]
then
FL="lista"
else
FL="$2"
fi
if [ "$3" == "" ]
then
SKIP=""
else
SKIP="-s $3"
fi
#$PERC/calcACV $FL 1000
$PERC/calcfqtself_linear $SKIP $FL 1000 $MAX $MAX 
$PERC/calcrho $FL $MAX $MAX
$PERC/calcfqtcoll $MAX $MAX 1000
