cat $1 | awk  '{ if ($3=="(SLOPE)") print $5}' > _allcoeff.tmp
MEAN=`cat _allcoeff.tmp | awk 'BEGIN { sum=0; i=0;} {sum=sum+$1; i=i+1;} END {print sum/i}'`
cat _allcoeff.tmp | awk -v mean=$MEAN 'BEGIN {i=0; sum=0;} {sum=sum+($1-mean)*($1-mean); i=i+1;} END {print ("mean=",mean,"stdv=",sqrt(sum/i),"abs error=",sqrt(sum/i)/sqrt(i),"rel err=",sqrt(sum/i)/sqrt(i)/mean);}'
