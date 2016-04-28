for f in uI_*
do
cd $f
cp ../start_NOCROWD_M* .
../create_dirs.sh
cd ..
done
