file=$1
cut=$2
#cut=0.1

cat $file |awk -v CUT=$cut '{if($2>CUT) print $1,$2}' > $file.xxx
#cat $file |awk -v CUT=$cut '{if($2<CUT) print $1,$2}' > $file.xxx

rm -f gnu.x
echo "g(x)=exp(-6.0*Drot*x)"  >> gnu.x
echo "fit g(x) '$file.xxx' via Drot" >> gnu.x
echo "beta = 1.01" >> gnu.x
echo "plateau=0.01"  >> gnu.x
echo "f(x) = plateau+exp (-(6*Drot*x)**beta)" >> gnu.x
echo "fit f(x) '$file.xxx' via Drot,beta,plateau" >> gnu.x
echo "print Drot,beta,plateau" >> gnu.x
gnuplot gnu.x >& aaa
tail -1 aaa; rm -f aaa
rm -f gnu.x $file.xxx
