#per usare questo script in mcsim_hardell va definita la macro TEST_DELAUNAY e va ricompilato
#inoltre vanno generate delle configurazioni che si devono trovare nella stessa directory dello script
DIR=$HOME
for f in `ls cnf-*`
do
  cat $f | gawk '{if (NR > 1) print ($1,$2)}' > $DIR/conf-python
  cp $f startcnf
  ./mcsim_hardell 2 > _aaa_
  cat _aaa_ | gawk '{if (NR >=7) {for (i=1; i <= NF; i++) {printf("%s",$i); if (i<NF) {printf(" ");} else {printf("\n");};};}}' > $DIR/nnmy.dat
  rm _aaa_
  python3 ./delaunay.py 
done
