# program to produce Mean Square Displacements
#
# calculates the MSD for all the files withe extension ".file.list"
# present in the directory. Such files are in the format:
#	number of configurations
#	directory of the configurations
#	first configuration
#	....
#	....
#	last configuration
#
# to produce these files, one can use the script makefilelist.sh
#
#

q=$1
l=$2
m=$3

make clean; 
make CORR-SQY;
echo GEOMETRIC MEAN 

ls -1 *.file.list|awk -F ".file.list" -v suff="S"$q"Y"$l$m '{print $1".file.list", $1""suff}'>abc
cat abc |while read list target
do
	./CORR-SQY.x $q $l $m < $list > $target
done
make clean
#rm abc

