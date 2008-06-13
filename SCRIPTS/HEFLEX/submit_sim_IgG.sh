if [ "$1" == "" ]
then
echo "Devi fornire la conf iniziale"
exit
fi
CONFDIR="$HOME/IgG/IgG_CONF"
EXEB="$HOME/IgG/MDsimul/bin/ellipsoid"
cd $HOME/IgG/ANTIBODYGENE
CDIR=`pwd`
DIR=`basename $CDIR`
PARFILE="ellipsoid_flex.par"
if [ "$2" == "" ]
then
TEMP="0.2"
SUFFIX=""
else 
SUFFIX="_T$2"
TEMP=$2
fi
NUM_ANTIBODY="200"
RESC_TIME="50.0"
ENESTPS="1000"
STORERATE="1000.0"
AGD=`echo $1 | awk -F - '{print $2}'`
ABD=`echo $1 | awk -F - '{print $3}'`
#RUN NUMBER
NR=`echo $1 | awk -F - '{print $4}'` 
if [ "$NR" != "" ]
then
NRSUFFIX="R"$NR 
fi
EEXE=`echo "ellips-"IgG"${ABD}-sig"${AGD}${NRSUFFIX}${SUFFIX}` 
ABDIR="IgG_$ABD"
if [ ! -e "$ABDIR" ] 
then
mkdir $ABDIR
fi
cd $ABDIR
AGDIR="sigma_${AGD}${NRSUFFIX}${SUFFIX}"
if [ ! -e "$AGDIR" ]
then
mkdir $AGDIR
fi
cp ../$PARFILE $AGDIR
cd $AGDIR
cp $CONFDIR/$1 .
cp $1 CorIni
rm -fr Store*
rm -fr COORD_TMP*
cp $EXEB .
ln -sf ./ellipsoid $EEXE 
../../set_params.sh $PARFILE stepnum 100000000 storerate $STORERATE stripStore 1 intervalSum 10.0 scalevel 1 rescaleTime $RESC_TIME hardwall 1 par2save 0-$[$NUM_ANTIBODY-1] inifile CorIni endfile corfinal.cor temperat $TEMP VSteps $ENESTPS
cat ../../template_pbs_submit_script.sh | awk -v LF="$CDIR/$ABDIR/$AGDIR" -v eexe="./$EEXE" -v jn="$EEXE" '{for (i=1; i <= NF; i++) {if ($i=="_LOCAL_FOLDER_") printf("%s/* ", LF); else if ($i=="_JOBNAME_") {printf("%s ", jn)} else if ($i=="_EXE_NAME_") {printf("%s ",eexe); } else printf("%s ", $i); } printf("\n");}' > pbs_submit_script.sh
#./$EEXE -fa ellipsoid_flex.par > screen &
#get last job ID
#JOBID=`jobs -p`
#write to a FILE
#echo $JOBID > ../../_JOBID_ 
