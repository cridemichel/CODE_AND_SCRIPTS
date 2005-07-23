grep "0 0 0 0" sq0* | awk '{print $6,$7}' > sq0000.dat
grep "2 0 2 0" sq0* | awk '{print $6,$7}' > sq2020.dat
grep "0 0 1 0" sq0* | awk '{print $6,$7,$8}' > sq1000.dat
grep "0 0 2 0" sq0* | awk '{print $6,$7,$8}' > sq2000.dat
grep "1-1 1-1" sq0* | awk '{print $6,$7,$8}' > sq1m11m1.dat
grep "1-1 2-1" sq0* | awk '{print $6,$7,$8}' > sq2m11m1.dat
grep "1 0 1 0" sq0* | awk '{print $6,$7,$8}' > sq1010.dat
grep "1 0 2 0" sq0* | awk '{print $6,$7,$8}' > sq2010.dat
grep "1 1 1 1" sq0* | awk '{print $6,$7,$8}' > sq1111.dat
grep "2 1 1 1" sq0* | awk '{print $6,$7,$8}' > sq2111.dat
grep "2-2 2-2" sq0* | awk '{print $6,$7,$8}' > sq2m22m2.dat
grep "2 1 2 1" sq0* | awk '{print $6,$7,$8}' > sq2121.dat
grep "2 2 2 2" sq0* | awk '{print $6,$7,$8}' > sq2222.dat
