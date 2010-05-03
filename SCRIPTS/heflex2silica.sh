NN=`cat $1 | awk -F : '{if ($1=="parnum") print $2}'`
#NA=`cat $1 | awk -F : '{if ($1=="parnumA") print $2}'`
NN="47125"
NA="250"
#echo "N=" $NN "NA= " $NA
cat $1 | gawk -v N=$NN -v NA="$NA" 'BEGIN {at=0; cc=0; first=1;} {if (at==2 && first) {print ("parnum:",N); print ("parnumA:", NA); print ("@@@"); first=0;} if ($0=="@@@") at++;}' 
cat $1 | gawk -v N=$NN -v NA="$NA" 'BEGIN {at=0; cc=0;} {if (at==2) { if (cc < N && $13=="0") print($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12); cc++; } if ($0=="@@@") at++;}'
cat $1 | gawk -v N=$NN -v NA="$NA" 'BEGIN {at=0; cc=0;} {if (at==2) { if (cc < N && $13=="1") print($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12); else if (cc >= N  && cc < 2*N) print($1,$2,$3,"0 0 0"); else if (cc==2*N) print $0; cc++; } if ($0=="@@@") at++;}'
