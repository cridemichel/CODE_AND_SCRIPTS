file=$1
cut=$2

cat $file |awk -v CUT=$cut '{if($2>CUT) print $1,$2}' > $file.xxx
rm -f gnu.x
echo "f(x) = exp (-6*Drot*x)" >> gnu.x
echo "fit f(x) '$file.xxx' via Drot" >> gnu.x
echo "plot f(x),'$file.xxx'" >> gnu.x
echo "pause -1" >> gnu.x
gnuplot gnu.x 
rm -f gnu.x $file.xxx
