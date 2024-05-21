for k in `ls -d X0_*`
do
 cd $k
 ../create_dirs.py
 cd ..
done
