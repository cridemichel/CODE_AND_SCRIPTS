if [ "$2" == "" ]
then
echo "please provide max and min N"
exit
else
NMIN="$1"
NMAX="$2"
fi
if [ "$3" == "" ]
then
NB="2"
else
NB="$1"
fi
TOTN="echo $NMAX-$NMIN+1| bc"
BL=`echo $TOTN/$NB| bc`
i=0
NI=0
while [ $i -le $NB ]  
do
NF=$[$NI+$BL]
if [ $NF -gt $NMAX ]
then
NF=$NMAX
fi
./restart_sus.sh $NI $NF 2>&1 > screen-$i &
NI=$[$N+$BL]
done
