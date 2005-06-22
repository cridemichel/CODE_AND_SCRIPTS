DIR=$1
#ls -1 $DIR |awk -v dir=$DIR '{print dir,$1}' > dirlist.dat
ls -1 $DIR |grep Phi |awk -v dir=$DIR '{print dir,$1}' > dirlist.dat


PWD=`pwd`
cat dirlist.dat | while read PREFIX DIR
do
	sh makefilelist-local.sh $PWD $PREFIX $DIR
done

