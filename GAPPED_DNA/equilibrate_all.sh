PF="ellipsoid_flex_mc.par"
for f in C*
do
cd $f
#cp start_phi* start.cnf
../set_params.sh $PF temperat 10.0 stepnum 100000
EXE=`ls gapdna*`
nohup mosrun ./$EXE -fa ellipsoid_flex_mc.par > screen &
sleep 1
cd ..
done
