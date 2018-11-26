MATDIMS="3 4 5 6 8 10 15 20 30 50 80 120 150 200" 
echo -n "" > results.out
echo -n "" > speedup.out
SF="test_lib.c++"
EXE="test_lib"
TC="gtime -f %e -o _tempo_ "
NCINI="2000000"
NC=$NCINI
DEFINES="-DTEST_INV "
CC="g++"
LIBS="-llapack -lblas -larmadillo"
FLAGS="-std=c++17 -ffast-math -O3"
for N in `echo $MATDIMS`
do
NC=`echo "e(2.2*l(3.0/${N}))*${NCINI}" | bc -l | LANG=C awk '{printf("%d", $0)}'`  
echo "Doing " $NC "loops for N=" $N
#my class
$CC -DNMAT="$N" -DSEL="1" $DEFINES $FLAGS $LIBS $SF -o $EXE
$TC ./$EXE $NC
TMY=`cat _tempo_`
#armadillo 
$CC -DNMAT="$N" -DSEL="3" $DEFINES $FLAGS $LIBS $SF -o $EXE 
$TC ./$EXE $NC 
TAR=`cat _tempo_`
#eigen
$CC -DNMAT="$N" -DSEL="4" $DEFINES $FLAGS $LIBS $SF -o $EXE
$TC ./$EXE $NC 
TEI=`cat _tempo_`
echo $N " " $TMY " " $TAR " " $TEI >> results.out
SUAR=`echo "${TAR}/${TMY}" | bc -l` 
SUEI=`echo "${TEI}/${TMY}" | bc -l` 
echo $N " " $SUAR " " $SUEI " " >> speedup.out  
done  
