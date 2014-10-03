if [ "$1" = "" ]
then
M=5
else
M=$1
fi
ls Qnc*  > lista
cat lista | sort -t - -k 4 -n | awk -v M="$M" 'BEGIN { b=0; i=0;} { if (b%M==0)\
 print $0;
 if (i==15)\
  {i=0;b=b+1;}\
 else\
  i=i+1;}'
