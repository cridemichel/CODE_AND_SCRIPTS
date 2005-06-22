file=$1
cut=$2
#cut=0.1

cat $file |awk -v CUT=$cut '{if($2>CUT) print $1,$2}' > $file.xxx
#cat $file |awk -v CUT=$cut '{if($2<CUT) print $1,$2}' > $file.xxx

rm -f gnu.x
echo "g(x)=exp(-6.0*Drot*x)"  >> gnu.x
echo "fit g(x) '$file.xxx' via Drot" >> gnu.x
echo "tau1 = 1.0/(6.0*Drot) "  >> gnu.x
echo "tau2 = 10*tau1 "  >> gnu.x
echo "beta = 2.0" >> gnu.x
echo "f(x) = cos(x/tau1)*(x/tau2)**-beta" >> gnu.x
echo "fit f(x) '$file.xxx' via tau1,tau2,beta" >> gnu.x
echo "plot f(x),'$file.xxx'" >> gnu.x
echo "pause -1" >> gnu.x
gnuplot gnu.x 
rm -f gnu.x $file.xxx
