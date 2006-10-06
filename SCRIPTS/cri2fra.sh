echo " Ricorda che serve il file ciclo.dat "
#
# da file Store a file Cnf
#
gunzip -q $1
echo $1 > junk
awk -F ".gz" '{print $1}' junk > junk2
nomefile=`cat junk2`
echo $nomefile
NPC=`awk '{print $2}' ciclo.dat`
echo "NPC="  $NPC
tempoconf=`awk -F "-" -v NPC=$NPC '{print $2*NPC+$3+1}' junk2 `
echo " Numero Configurazione="  $tempoconf
BX=`tail -1 $nomefile`
echo " BX= " $BX
time=` awk -F ":" '{if ($1=="time") printf " %i \n ", int($2*10000)}'  $nomefile `
echo " TIME=" $time
echo  $time $time 350 350 $tempoconf > head
#
#  il quinto numero qui dovrebbe essere il Dt -- bilancia il 100000 usato sopra
#
echo $BX $BX $BX  0 0.0001 0 >> head
#
# qui ricostruisce la molecola come H H O (migliorare il printf)
#
awk -F":" 'BEGIN {cc=0} {if (cc==2) print $0; if ($1=="@@@") cc=cc+1; }' $nomefile > yyy
 awk  '{if (NR<351) printf "%f %f %f \n %f %f %f \n %f %f %f \n " , $1+$4, $2+$5, $3+$6, $1+$7, $2+$8, $3+$9, $1, $2, $3 }' yyy > xxx
cat head xxx >  Cnfpwm-$tempoconf


