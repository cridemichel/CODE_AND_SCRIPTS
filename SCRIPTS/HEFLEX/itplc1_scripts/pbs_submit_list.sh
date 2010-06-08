for f in `cat $1`
do
./submit_sim_IgG.sh $f 0.2 $2
done
