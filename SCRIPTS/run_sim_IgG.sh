if [ "$1" == "" ]
then
echo "Devi fornire la conf iniziale"
exit
fi
CONFDIR="$HOME/IgG_CONF"
CDIR=`pwd`
DIR=`basename $CDIR`
EEXE=`echo "ellips"$DIR$1` 
PARFILE="ellipsoid_flex.par"
if [ $2 == "" ]
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
STORERATE="100.0"
AGD=`echo $1 | awk -F - '{print $2}'`
ABD=`echo $1 | awk -F - '{print $3}'`
ABDIR="IgG_$ABD"
if [ ! -e "$ABDIR" ] 
then
mkdir $ABDIR
fi
cd $ABDIR
AGDIR="sigma_$AGD$SUFFIX"
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
ln -sf ../../ellipsoid $EEXE 
../../set_params.sh $PARFILE stepnum 100000000 storerate $STORERATE stripStore 1 intervalSum 10.0 scalevel 1 rescaleTime $RESC_TIME hardwall 1 par2save 0-$[$NUM_ANTIBODY-1] inifile CorIni endfile corfinal.cor temperat $TEMP VSteps $ENESTPS
$EEXE -fa ellipsoid_flex.par > screen &
