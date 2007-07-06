echo "Job started..."
FC=0
ls Store* | sort -t - -k 2 -k 3 -n > _sorted.tmp
for f in `cat _sorted.tmp`
do 
$HOME/molgl/molgl -di -70.0 -sq -ni -f frame.${FC}.png $f
FC=$[$FC+1]
done
rm _sorted.tmp
ls frame.*.png | sort -t '.' -k 2 -n > _sorted.tmp
convert `cat _sorted.tmp` protmovie.mpeg
echo "...done!"
