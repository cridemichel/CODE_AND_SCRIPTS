for f in X0_*
do
cd $f
echo "calculating" $f
cat allTAUtra-${f}.dat | awk '{print $2,$3}' >  tautra.dat
A=`cat tautra.dat | awk '{print $1}'|sort -k 1 -n | head -1`
B=`cat tautra.dat | awk '{print $1}'|sort -k 1 -n | tail -1`
echo "A="$A " B=" $B
echo "p1=" $A "; p2=" $B "; fit [p1:p2] a*x**3+b*x**2+c*x+d \"tautra.dat\" via a,b,c,d; print p1, p2, a, b, c, d" > fit.tmp
echo "Creo la spline di tautra.dat"
#gnuplot fit.tmp > gpout.tmp 2>&1  
#tail -1 gpout.tmp > cubic.par
#../gencubicfit cubic.par > tautra_fit.dat
echo "READ \"tautra.dat\"" > spline.agrs
echo "INTERPOLATE(S0,MESH(0.15,0.60,100),ASPLINE,ON)" >> spline.agrs
echo "WRITE G0.S1 FILE \"tautra_spline.dat\"" >> spline.agrs
echo "EXIT" >> spline.agrs
xmgrace  -noask -nosafe -batch spline.agrs
cd ..
done
