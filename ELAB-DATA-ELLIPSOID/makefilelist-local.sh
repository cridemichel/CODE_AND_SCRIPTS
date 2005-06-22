# use in a directory; it produces in ../ a file with the name 
# (name of the directory).file.list in such a format:
#
#       number of configurations
#       directory of the configurations
#       first configuration
#       ....
#       ....
#       last configuration
#
# those files are normally used to calculate averages over the  
# configurations
#

MYPWD=$1
PREFIX=$2
DIR=$3

cd $PREFIX/$DIR
echo ---; pwd; echo ---;
echo writing to the file $MYPWD/$DIR.file.list; echo;echo;

gunzip Store-*-0.gz

ls -1 Store-*-0 | wc -l > $MYPWD/$DIR.file.list
echo $PREFIX/$DIR/ >> $MYPWD/$DIR.file.list
ls -1  Store-?-0    >> $MYPWD/$DIR.file.list
ls -1  Store-??-0   >> $MYPWD/$DIR.file.list
ls -1  Store-???-0  >> $MYPWD/$DIR.file.list
ls -1  Store-????-0 >> $MYPWD/$DIR.file.list

cd $MYPWD
