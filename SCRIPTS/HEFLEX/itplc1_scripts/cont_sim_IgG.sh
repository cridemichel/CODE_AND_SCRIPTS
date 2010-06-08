if [ "$1" == "" ]
then
echo "Devi fornire la directory iniziale della simulazione da riavviare e il nome dell'eseguibile"
exit
fi
CDIR=`pwd`
cd $1
cat ../../template_pbs_cont_script.sh | awk -v LF="$CDIR/$1" -v eexe="./$2" -v jn="$2" '{for (i=1; i <= NF; i++) {if ($i=="_LOCAL_FOLDER_") printf("%s/ ", LF); else if ($i=="_JOBNAME_") {printf("%s ", jn)} else if ($i=="_EXE_NAME_") {printf("%s ",eexe); } else printf("%s ", $i); } printf("\n");}' > pbs_cont_script.sh
rm *.e* > /dev/null 2>&1
rm *.o* > /dev/null 2>&1
qsub -q low pbs_cont_script.sh
