#!/bin/sh
PROGNAME="basins"                 # nome dell'eseguibile da riavviare
NUM_OF_NODES=3                    # numero di nodi del run parallelo (n. pc in ~/boot)
LAMFILE=${HOME}/boot              # file per lamboot
EXEPERC=${HOME}/MDsimul/bin     # percorso del file eseguibile da riavviare
PERC=${HOME}/analisi_dati/cmprmk80/M14PRn # percorso delle config. e del file files.dat
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
PATH=$PATH:/usr/local/lam-6.3.2/bin/:${HOME}:.
LAMHOME=/usr/local/lam-6.3.2
# ==================================================================================
export PATH LAMHOME
cd $PERC
ps axl > processi.a
# ora il test dell'esistenza del programma eseguiti con lamboot è piu' accurato, infatti
# verifica che esistano tutti i processi ("mpirun..basins", "lamd..", 2 x "basins") 
# dell'utente che sta lanciando riavvia.sh
# echo "UID:" $UID
if [ "$1" = "lam" -o "$1" = "" ]
then
PROG_EXIST=`awk -v PN=$PROGNAME -v UID=$UID 'BEGIN {count = 0}; $2==UID && $0 ~ PN && /mpirun/ {count++}; $2==UID && /lamd/ && ! /mpirun/ {count++}; $2==UID && $0 ~ PN && ! /mpirun/ {count++}; END { if (count < 4) {print 0} else {print count} } ' processi.a`
else
# nel caso del programma seriale ci sara' solo un'occorrenza dello stesso...
PROG_EXIST=`awk -v PN=$PROGNAME -v UID=$UID 'BEGIN {count = 0}; $2==UID && /PN/ {count++}; END { if (count < 1) {print 0} else {print count} } ' processi.a`
fi
#echo $PROG_EXIST
#exit
rm processi.a
ls CnfT*-* > cnf.a	
ls QncT*-* 2>/dev/null 1> qnc.a
if [ $[PROG_EXIST] -eq 0 ]
then
    #il programma non sta girando dunque va rilanciato 
    echo "Riavvio in corso..."
    mancanti cnf.a qnc.a > todo.a
    if [ `cat todo.a | wc -c` -le 1 ]
    then
	echo "Finito!"
	exit
    fi
    echo `pwd`/ > files.dat ; cat todo.a >> files.dat
    if [ "$1" = "lam" -o "$1" = "" ]
    then
	${LAMHOME}/bin/wipe -v ${LAMFILE} # uccide tutti processi lam in giro
	lamboot -v ${LAMFILE} # rilancia lamd
	mpirun -lamd -c $NUM_OF_NODES -s h ${EXEPERC}/$PROGNAME -- ${PERC}/files.dat  &
    else
	${EXEPERC}/$PROGNAME ${PERC}/files.dat &
    fi
else
    if [ "$1" = "lam" -o "$1" = "" ]
    then
	# stampa anche lo stato di tutti i task usando mpitask 
	${LAMHOME}/bin/mpitask
    fi
    echo "------------------------------------------------------------------------------- "
    echo "Working Directory: " ${PERC}/
    echo `mancanti cnf.a qnc.a >/dev/null`
fi





