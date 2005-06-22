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

make clean; 
make CORR-Y;
echo GEOMETRIC MEAN 

ls -1 *.file.list|awk -F ".file.list" '{print $1".file.list", $1".corrY"}'>abc
cat abc |while read list target
do
	./CORR-Y.x < $list > $target
done
make clean
rm abc

