file=$1
rm -f gnu.x
echo "f(x) = 6*D*x + r0" >> gnu.x
echo "fit f(x) '$file' via D,r0" >> gnu.x
echo "plot f(x),'$file'" >> gnu.x
echo "pause -1" >> gnu.x
gnuplot gnu.x 
rm -f gnu.x
