EL=$1
PHI=$2

DIR="../nodo14/scala/GRID_N256/SAVE-N256-EL$EL/Phi"$PHI"_"$EL

gunzip $DIR/Store-*.gz

echo `ls -1 $DIR/Store-* |wc -l` > aaa
echo $DIR/ >> aaa
ls -1 $DIR/ |grep Store >> aaa

make Y20Y20
./Y20Y20.x < aaa

gzip $DIR/Store-*

