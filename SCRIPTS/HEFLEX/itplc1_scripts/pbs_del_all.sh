JL=`qstat -a | grep ellips | awk '{print $1}' | awk -F '.' '{print $1}'`
for f in $JL
do
echo "Killing: " $f
qdel $f
done
