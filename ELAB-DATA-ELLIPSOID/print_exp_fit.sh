file=$1
cut=$2
#cut=0.1

cat $file |awk -v CUT=$cut '{if($2>CUT) print $1,$2}' > $file.xxx
rm -f gnu.x
echo "f(x) = exp (-6*Drot*x)" >> gnu.x
echo "fit f(x) '$file.xxx' via Drot" >> gnu.x
echo "print Drot" >> gnu.x
gnuplot gnu.x >& aaa
tail -1 aaa; rm -f aaa
rm -f gnu.x $file.xxx
