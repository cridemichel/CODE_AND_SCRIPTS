MATDIMS="4 5 6 7" 
#MATDIMS=" 3 4 5 6 8 10 15 20 30 50 80 120 150 200"
echo -n "" > results.out
echo -n "" > speedup.out
SF="test_lib.c++"
EXE="test_lib"
#TC="gtime -f %e -o _tempo_ "
TC="time -p "
NCINI="1000000"
NC=$NCINI
DEFINES="-DTEST_MUL "
ENABLE_GLM="0"
if [ "$1" == "" ]
then
TESTS="all"
else
TESTS="$1"
fi
CC="g++ -framework accelerate"
#CC="g++-8  -I/usr/local/include -L/usr/local/lib -L/usr/local/opt/openblas/lib/ -L/usr/local/opt/lapack/lib/"
LIBS="-llapack -lblas -larmadillo"
FLAGS="-std=c++17 -march=native -ffast-math -O3"
#FLAGS="-std=c++17 -ffast-math -O3 "
for N in `echo $MATDIMS`
do
NC=`echo "e(2.2*l(3.0/${N}))*${NCINI}" | bc -l | LANG=C awk '{printf("%d", $0)}'`  
echo "Doing " $NC "loops for N=" $N
#my class
if [ \( "$TESTS" == "all" \) -o \( "$TESTS" == "MY" \) ]
then
$CC -DNMAT="$N" -DSEL="1" $DEFINES $FLAGS $LIBS -o $EXE $SF 
$TC ./$EXE $NC 2>&1 | LANG=C awk '/real/{print $2}' > _tempo_
TMY=`cat _tempo_`
fi
#glm
if [ \( \( "$TESTS" == "all" \) -o \( "$TESTS" == "GLM" \) \) -a \( "$ENABLE_GLM" == "1" \) ]
then
$CC -DNMAT="$N" -DSEL="2" $DEFINES $FLAGS $LIBS $SF -o $EXE 
$TC ./$EXE $NC 2>&1 | LANG=C awk '/real/{print $2}' > _tempo_
TGL=`cat _tempo_`
fi
#armadillo 
if [ \( "$TESTS" == "all" \) -o \( "$TESTS" == "ARM" \) ]
then
$CC -DNMAT="$N" -DSEL="3" $DEFINES $FLAGS $LIBS $SF -o $EXE 
$TC ./$EXE $NC 2>&1 | LANG=C awk '/real/{print $2}' > _tempo_
TAR=`cat _tempo_`
fi
#eigen
if [ \( "$TESTS" == "all" \) -o \( "$TESTS" == "EIG" \) ]
then
$CC -DNMAT="$N" -DSEL="4" $DEFINES $FLAGS $LIBS $SF -o $EXE
$TC ./$EXE $NC 2>&1 | LANG=C awk '/real/{print $2}' > _tempo_
TEI=`cat _tempo_`
fi
if [ "$TESTS" == "all" ]
then
if [ "$ENABLE_GLM" == "1" ]
then
echo $N " " $TMY " " $TGL " " $TAR " " $TEI >> results.out
else
echo $N " " $TMY " " $TAR " " $TEI >> results.out
fi
SUAR=`echo "${TAR}/${TMY}" | bc -l` 
SUEI=`echo "${TEI}/${TMY}" | bc -l` 
if [ "$ENABLE_GLM" == "1" ]
then
SUGL=`echo "${TGL}/${TMY}" | bc -l` 
echo $N " " $SUGL " " $SUAR " " $SUEI " " >> speedup.out  
else
echo $N " " $SUAR " " $SUEI " " >> speedup.out  
fi
fi
if [ \( "$TESTS" == "GLM" \) -a \( "$ENABLE_GLM" == "1" \) ]
then
echo $N " " $TGL >> results.out
fi
if [ "$TESTS" == "MY" ]
then
echo $N " " $TMY >> results.out
fi
if [ "$TESTS" == "ARM" ]
then
echo $N " " $TAR >> results.out
fi
if [ "$TESTS" == "EIG" ]
then
echo $N " " $TEI >> results.out
fi
done
