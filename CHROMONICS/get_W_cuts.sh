cat Welconst.dat | LANG=C awk '{if ($1=="5.000000") print ($2,$3,$4)}' > W13_W23_vs_k3.dat
cat W13_W23_vs_k3.dat | LANG=C awk '{print ($1,$2)}' > W13_vs_k3.dat
cat W13_W23_vs_k3.dat | LANG=C awk '{print ($1,$3)}' > W23_vs_k3.dat
cat Welconst.dat | LANG=C awk '{if ($2=="5.000000") print ($1,$3,$4)}' > W13_W23_vs_k1.dat
cat W13_W23_vs_k1.dat | LANG=C awk '{print ($1,$2)}' > W13_vs_k1.dat
cat W13_W23_vs_k1.dat | LANG=C awk '{print ($1,$3)}' > W23_vs_k1.dat
