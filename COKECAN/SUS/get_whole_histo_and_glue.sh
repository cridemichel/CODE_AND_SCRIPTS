if [ "$1" == "" ]
then
ls -d N_*_* | sort -t _ -k 2 -n > lista
./glue_histo lista
else
./glue_histo $1
fi
#echo -n "" > histo.dat
#for f in `cat lista`
#do 
#cat $f/histo.dat >> histo.dat
#done
