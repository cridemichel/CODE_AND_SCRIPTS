ls -1 STORE-2.0/Store-*-* | head -1 |awk -F "/" '{print $1,$2}' > aaa.in

make clean; make SHOW

cc interval.c -o interval.x

./SHOW.x < aaa.in |awk '{print $3,$4}' | ./interval.x > qmag.list
 
rm -f aaa.in

make clean
