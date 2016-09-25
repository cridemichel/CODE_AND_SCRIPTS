#!/bin/sh
MEASDIR=$HOME/simdat/measures
BASENAME="BMsoftT50_PP8NVE"
MDEXE=$HOME/MDsimul/bin/bimix
function md_subst_one()
{
	#questa funzione sostituisce il valore del parametro $2 con il valore $3 nel file $1
	awk -F : -v parval=$3 -v parlbl=$2  '$1==parlbl {print (parlbl ": " parval)} $1 != parlbl { print $0 }' $1 > /tmp/md_subst.tmp
	cp /tmp/md_subst.tmp $1 
	rm /tmp/md_subst.tmp
}
if [ "$1" == "" ]
then
	echo "No volumes file supplied!"
exit 1
fi
function md_execall()
{
#La prima linea deve contenere la directory dove si trovano i .par files
#La seconda linea deve contenere i nomi dei file di parametri
#La terza linea deve contenere i nomi dei parametri
#inoltre la prima colonna è riservata al numero che indica il file di parametri
#da usare e la seconda colonna serve per indicare la directory in cui salvare le misure
read PARDIR
read -a PARFILES
i=0
while [ "${PARFILES[$i]}" != "" ] 
do 
	i=$((i+1))
done
NUM_PF=$i
if [ $((i)) -eq 0 ] 
then
	echo "Devi specificare almeno un file di parametri..."
	exit 1
fi
read -a PARS
#le linee successive contengono i valori con cui fare l'equilibratura e la produzione
read -a PARSVAL
while [ "${PARSVAL[0]}" != "" ] 
do
	if [ ${PARSVAL[0]} -gt $((NUM_PF)) -o ${PARSVAL[0]} -lt 0 ]
	then
		echo "Il file di paramentri n."  $NUM_PF " non è stato specificato!"
		exit 1
	fi
	if [ ${PARSVAL[0]} -gt 0 ]
	then
		PFILE=${PARDIR}/${PARFILES[$((${PARSVAL[0]}-1))]}
		echo "Using parfile: " $PFILE
		if [ ${PARSVAL[1]} == "auto" ] 
		then
			BN="${MEASDIR}/$BASENAME"
			MDIRYN="auto"
		elif [  ${PARSVAL[1]} == "none" ]
		then	
			MDIRYN="none"
		else
			BN=${MEASDIR}/${PARSVAL[1]}
			MDIRYN="yes"
		fi
		i=2
		while [ "${PARS[$i]}" != "" ]
		do
			echo "Setting par: " ${PARS[$i]}
			md_subst_one $PFILE ${PARS[$i]} ${PARSVAL[$i]}
			if [ "$MDIRYN" == "auto" ]
			then
				BN="${BN}_${PARS[$i]}${PARSVAL[$i]}" 
			fi
			i=$((i+1))
		done
		#Running
		$1 -f  $PFILE
		if [ $MDIRYN != "none" ]
		then
			if [ ! -d "$BN" ] 
			then
				mkdir $BN
			else
				rm -f $BN/* 
			fi
		cp -f ${MEASDIR}/*-0 $BN
		fi
	fi
	read -a PARSVAL
done
}
#NOTA:
#il file con i set di parametri deve essere del tipo:
# <percorso parfiles>
# <file1.par>	<file2.par> ....
# ParFile 	MeasDir		<parameter label #1>	<parameter label #2> ...
# <num>		<measdir>	value			value
# .
# .
# .
# dove <num> è un numero intero che indica il file di parametri da usare per la 
# simulazione corrente e <measdir> = "auto" per far generare automaticamente il nome 
# della directory oppure "none" per non salvare i file di misura, oppure è 
# la directory in cui salvare i file di misura (senza specificare il percorso!) 
#Execute simulation
cat $1 | md_execall $MDEXE 
