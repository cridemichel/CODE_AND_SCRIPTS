for f in C*
do
cd $f
cp start.cnf start_bak.cnf
CF=`ls -rt COORD_TMP_ASCII*| tail -1`
cat $CF | awk 'BEGIN{at=0} {if (at >= 1) print $0; if ($1=="@@@") at++;}' > start.cnf
cd ..
done
