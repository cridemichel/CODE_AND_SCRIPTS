PR="$HOME/HSMIX/SIMULATIONS/SCRIPTS/calcScatFunc_allq.sh"
for f in `cat $1`
do
cd $f
ls Store-*-*| sort -t - -k 2 -k 3 -n > lista
$PR &
cd ..
done


