cat Welconst.dat  | LANG=C gawk '{if ($1*$1 < 5.5 && $2*$2 < 5.5) printf("%.15G %.15G %.15G\n",$1*$1,$2*$2,$4)}' > W23.dat

