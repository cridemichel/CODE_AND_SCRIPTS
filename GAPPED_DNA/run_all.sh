PF="ellipsoid_flex_mc.par"
for f in C*
do
cd $f
#cp start_phi* start.cnf
#cp CorFinal start.cnf
../set_params.sh $PF temperat 0.12 stepnum 800000000 bakStepsAscii 20000
EXE=`ls gapdna*`
nohup mosrun ./$EXE -fa ellipsoid_flex_mc.par > screen &
sleep 1
cd ..
done
