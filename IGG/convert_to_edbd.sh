for f in `cat $1`
do
cat $f | awk 'BEGIN{IGN=-1; at=0;} {if ($0=="@@@") at++; if ($1=="250" && $2="250" && $3="250" && $4="250" && $5="1000") print ("250 250 250 250",$6); else if ($1=="ntypes:") print("ntypes: 5"); else if ($1=="parnum:") print("parnum:",$2-1000); else if ($0=="1.000000 1.000000 1.000000 1.000000 0 1 ") print ("1.000000 1.000000 1.000000 1.000000 1 1"); else if ($0=="1.000000 1.000000 1.000000 1.000000 0 0 ") print ("1.000000 1.000000 1.000000 1.000000 1 0"); else if ($0=="0 1 5 0 1 0 1 1") print("0 1 4 0 1 0 1 1"); else if ($0=="1 1 5 0 1 0 1 1") print ("1 1 4 0 1 0 1 1"); else if ($1=="0.250" && $2=="0.250" && $3=="0.250") IGN=NR; else if (IGN==-1 ||NR - IGN > 3) print $0; if (at==2) exit;}' > conf_template
cat $f | awk 'BEGIN {beg=-1; lastnf=-1; at=0;} {nf=NF; if (at==2) {if (beg!=-1) {if (NR-beg < 1000 || NR-beg >= 2000) print $0; } else if ($13!="4" && $13!="5") print $0; else if ($13=="5") {for (n=1; n <=12; n++) printf("%s ",$n); printf("%d\n", 4)}}; if (nf==6 && lastnf==13){ beg=NR; };  lastnf=nf; if ($0=="@@@") at++;}' > conf_coords
cat conf_template  conf_coords > flexEDBD$f

done
rm conf_template conf_coords
