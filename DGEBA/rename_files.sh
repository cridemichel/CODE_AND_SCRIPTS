for f in `ls -d RUN_*`
do
cd $f
mv start.cnf store_file
cd ..
done
