#parameters for the spherical harmonics
L=$1
M=$2

#elongation to consider
EL=$3

#directories with the files
DIR=SQY$L.$M-STORE-$EL

echo "Doing averages in "$DIR

ls -1 $DIR/SQY* | wc  |awk '{print $1}' >a
`ls -1 $DIR/SQY* | head -1 | awk '{print "wc",$1}'` |awk '{print $1}' >b

paste a b > files$L.$M.list; rm -f a b

ls -1 $DIR/SQY* >> files$L.$M.list

rm -f AVE-SQY.x; make AVE-SQY

cat files$L.$M.list | ./AVE-SQY.x  > Sq-$EL-Y$L.$M.dat

rm -f AVE-SQY.x files$L.$M.list
