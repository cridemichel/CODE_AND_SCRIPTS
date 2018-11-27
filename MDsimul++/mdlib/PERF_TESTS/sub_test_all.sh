MATDIMS="2" 
echo -n "" > results.out
echo -n "" > speedup.out
SF="test_lib.c++"
EXE="test_lib"
TC="gtime -f %e -o _tempo_ "
NCINI="5000000"
NC=$NCINI
DEFINES="-DTEST_INV "
ENABLE_GLM="1"
if [ "$1" == "" ]
then
TESTS="all"
else
TESTS="$1"
fi
CC="g++"
LIBS="-llapack -lblas -larmadillo"
#FLAGS="-mfpmath=sse -std=c++17 -mavx512vbmi -ffast-math -O3"
FLAGS=" -std=c++17 -ffast-math -O3 "
for N in `echo $MATDIMS`
do
NC=`echo "e(2.2*l(3.0/${N}))*${NCINI}" | bc -l | LANG=C awk '{printf("%d", $0)}'`  
echo "Doing " $NC "loops for N=" $N
#my class
if [ \( "$TESTS" == "all" \) -o \( "$TESTS" == "MY" \) ]
then
$CC -DNMAT="$N" -DSEL="1" $DEFINES $FLAGS $LIBS -o $EXE $SF 
$TC ./$EXE $NC
TMY=`cat _tempo_`
fi
#glm
if [ \( \( "$TESTS" == "all" \) -o \( "$TESTS" == "GLM" \) \) -a \( "$ENABLE_GLM" == "1" \) ]
then
$CC -DNMAT="$N" -DSEL="2" $DEFINES $FLAGS $LIBS $SF -o $EXE 
$TC ./$EXE $NC 
TGL=`cat _tempo_`
fi
#armadillo 
if [ \( "$TESTS" == "all" \) -o \( "$TESTS" == "ARM" \) ]
then
$CC -DNMAT="$N" -DSEL="3" $DEFINES $FLAGS $LIBS $SF -o $EXE 
$TC ./$EXE $NC 
TAR=`cat _tempo_`
fi
#eigen
if [ \( "$TESTS" == "all" \) -o \( "$TESTS" == "EIG" \) ]
then
$CC -DNMAT="$N" -DSEL="4" $DEFINES $FLAGS $LIBS $SF -o $EXE
$TC ./$EXE $NC 
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
