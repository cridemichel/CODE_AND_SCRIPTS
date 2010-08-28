if [ "$3" == "" ]
then
MODE="loglog"
else
MODE="$3"
fi
if [ "$2" == "" ]
then
UGFACT="0.6"
else
UGFACT="$2"
fi
if [ "$UGFACT" == "-1.0" ]
then
UG="1.0"
else  
UG=`echo "${UGFACT}*0.129^2" | bc -l`
fi
#per awk log(x) = logaritmo naturale di x
if [ "$MODE" == "loglog" ]
then
cat $1 | LC_NUMERIC=C awk -v ug2=${UG} '{print (log(1.0/$1)/log(10.0)+log(ug2)/log(10.0),log($2)/log(10))}'
elif [ "$MODE" == "linlog" ]
then
cat $1 | LC_NUMERIC=C awk -v ug2=${UG} '{print (ug2/$1,log($2)/log(10))}'
elif [ "$MODE" == "linlin" ]
then
cat $1 | LC_NUMERIC=C awk -v ug2=${UG} '{print (ug2/$1,$2)}'
fi
