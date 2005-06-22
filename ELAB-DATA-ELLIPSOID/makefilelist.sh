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


DIR=`pwd`
ls -1 Store-*-0 | wc -l > $DIR.file.list
echo $DIR/ >> $DIR.file.list
ls -1  Store-?-0    >> $DIR.file.list
ls -1  Store-??-0   >> $DIR.file.list
ls -1  Store-???-0  >> $DIR.file.list
ls -1  Store-????-0 >> $DIR.file.list


