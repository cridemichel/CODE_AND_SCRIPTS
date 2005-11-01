gunzip Store*gz
ls Store* > junk
cat junk | while read nomefile
do
echo  $nomefile
REFTIME=`awk -F":" 'BEGIN {found=0} {if ($1=="refTime") { print $2; found=1}} END {if (!found) printf("0.0");}' $nomefile`
echo "reftime:" $REFTIME
awk -v reft="$REFTIME" -F":" 'BEGIN {cc=0} {if (cc>0) {if ($1=="time") printf("time: %.15G\n", reft+$2); else print $0}; if ($1=="@@@") cc=cc+1; }' $nomefile > Cnf$nomefile
#awk -F":" 'BEGIN {cc=0} {if (cc>0) print $0; if ($1=="@@@") cc=cc+1; }' $nomefile > Cnf$nomefile
done
gzip Store*
