echo "PRINT TO \"output.eps\"" > _aaa_
echo "DEVICE \"EPS\" OP \"level2\"" >> _aaa_
echo "PRINT" >> _aaa_
for f in $@
do
echo -n "processing file " $f "..." 
#LC_NUMERIC=C xmgrace -hdevice EPS -hardcopy $f
LC_NUMERIC=C xmgrace  -batch _aaa_ -hardcopy -nosafe $f
echo "done"
done
rm _aaa_
