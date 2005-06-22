cat phi_and_X0.dat |while read PHI EL
do
  ls -1 ABC/Phi"$PHI"_"$EL".msd | awk '{print "sh print_lin_fit.sh",$1}' > fst1
  sh fst1 >> fst2
done

echo "# phi X_0 D r0 (linear fit)" > msd.dat
paste phi_and_X0.dat fst2 >> msd.dat


rm fst?



