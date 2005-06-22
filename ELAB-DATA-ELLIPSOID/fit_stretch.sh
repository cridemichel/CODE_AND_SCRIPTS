
cat phi_and_X0.dat |while read PHI EL
do
  ls -1 ABC/Phi"$PHI"_"$EL".corrY | awk '{print "sh print_stretch_fit.sh",$1}' > fst1
  sh fst1 >> fst2
done

echo "# phi X_0 Drot beta plateau (stretched exponential fit)" > stretch.dat
paste phi_and_X0.dat fst2 >> stretch.dat

#cat fst3 | awk '{ if(($3<2)&&($3>0))   print $2,$3}' > tau_stretch.dat  
#cat fst3 | awk '{ if(($4<1.1)&&($4>0)) print $2,$4}' > beta_stretch.dat

#cat fst3 | awk '{ if(($3<2)&&($3>0))    print $1,$3}' > tau_stretch.dat  
#cat fst3 | awk '{ if(($4<1.1)&&($4>0.)) print $1,$4}' > beta_stretch.dat

#cat fst3 | awk '{ if(($3<2)&&($3>0))    print $2,$3,$1}' > tau_stretch.dat  
#cat fst3 | awk '{ if(($4<1.1)&&($4>0.)) print $2,$4,$1}' > beta_stretch.dat

rm fst?



