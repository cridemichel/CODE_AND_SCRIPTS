for f in `ls -1 -d sigE*`
do
cd $f
./calc_rates.sh &
sleep 0.5
cd ..
done
