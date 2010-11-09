echo "num particelle=" `cat $1 | LC_NUMERIC=C awk '{sum+=$1*$2;} END{print sum}'`
cat $1 | LC_NUMERIC=C awk '{sum+=$1*$2; sum2+=$2} END{print sum/sum2}'

