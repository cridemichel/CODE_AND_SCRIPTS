ls -d N_*_* | sort -t _ -k 2 -n > lista
#echo -n "" > histo.dat
#for f in `cat lista`
#do 
#cat $f/histo.dat >> histo.dat
#done
./glue_histo lista
