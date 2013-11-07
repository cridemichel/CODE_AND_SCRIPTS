echo -n "" > cov_vs_s.dat
ls -d *-* | sort -t - -k 1 -n > _lista_
for f in `cat _lista_`
do
echo "f=$f"
cd $f
COV=`cat covolume.dat | tail -1 | awk '{print $2}'`
S=`echo $f| awk -F - '{print $1}'`
echo $S $COV >> ../cov_vs_s.dat
cd ..
done
rm _lista_
