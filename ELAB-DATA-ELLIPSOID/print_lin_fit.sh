file=$1
rm -f gnu.x
echo "f(x) = 6*D*x + r0" >> gnu.x
echo "fit f(x) '$file' via D,r0" >> gnu.x
echo "print D,r0" >> gnu.x
gnuplot gnu.x >& aaa
tail -1 aaa; rm -f aaa
rm -f gnu.x
