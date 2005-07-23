ls Store* > junk
cat junk | while read nomefile
do
echo  $nomefile
awk -F":" 'BEGIN {cc=0} {if (cc>0) print $0; if ($1=="@@@") cc=cc+1; }' $nomefile > Cnf$nomefile
done

