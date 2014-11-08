for f in `cat $1`
do
PARNUM=`cat $f | awk  '{if ($1=="parnum:") print $2;}'`
cat $f | awk 'BEGIN{IGN=-1; at=0;} {if ($0=="@@@") at++; if ($1=="250" && $2=="250" && $3=="250" && $4=="250") print ("50 50 50 50",$5); else if ($1=="parnum:") print("parnum:",$2-800); else print $0; if (at==2) exit;}' > conf_template
#echo "PARNUM=" $PARNUM 
cat $f | awk -v NP=$PARNUM 'BEGIN {beg=-1; lastnf=-1; at=0;} {nf=NF; if (at==2 && beg==-1) beg=NR; if (at==2) {if ((NR-beg >=0 &&NR-beg < 200) || (NR-beg >= 1000 && NR-beg < NP) || (NR-beg-NP >= 0 && NR-beg-NP < 200) || (NR-beg-NP >= 1000 && NR-beg-NP < NP) ) print $0;}; if ($0=="@@@") {at++;}}' > conf_coords
tail -1 $f >> conf_coords
cat conf_template  conf_coords > RED-$f
done
rm conf_template conf_coords
