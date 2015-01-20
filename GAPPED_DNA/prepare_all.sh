PN="1000"
CONCS="150 210 230 250 271 285 300 320"
FACT="0.0008021977075921017"
for f in $CONCS
do
PHI=`echo $f*$FACT|bc -l|LANG=C gawk '{printf("%.4f\n",$1)}'`
echo "PHI=" $PHI " CONC=" $f " PARNUM=" $PN
./create_GDNA_conf start_phi${PHI}_c${f}.cnf $PN $PHI 2
done
