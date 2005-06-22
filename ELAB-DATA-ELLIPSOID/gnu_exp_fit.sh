file=$1
rm -f gnu.x
#echo "f(x) = mag * exp (-x/tau)" >> gnu.x
#echo "fit f(x) '$file' via mag,tau" >> gnu.x
echo "f(x) = exp (-x/tau)" >> gnu.x
echo "fit f(x) '$file' via tau" >> gnu.x
echo "plot f(x) , '$file'" >> gnu.x
echo "pause -1 \"Hit return to continue\"" >> gnu.x
echo "print tau" >> gnu.x
gnuplot gnu.x
rm -f gnu.x
