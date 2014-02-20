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
NB="$3"
fi
TOTN=`echo $NMAX-$NMIN+1 | bc`
BL=`echo $TOTN/$NB| bc`
echo "num blocks=" $NB " block length = " $BL
i=1
NI=$NMIN
FINE="0"
NF=$[$NI+$BL]
if [ $NI -ge $NMAX ]
then
echo "nothing to do, exiting..."
exit
fi
while [ $FINE == "0"  ]  
do
if [ $NF -gt $NMAX ]
then
NF=$NMAX
fi
if [ $i -eq 1 ]
then
WF="1"
else
WF="0"
fi
echo "block #" $i " NI=" $NI " NF=" $NF 
./restart_sus.sh $NI $NF $WF > screen-$i 2>&1 &
NI=$[$NI+$BL]
NF=$[$NI+$BL]
i=$[$[i]+1]
if [ $NI -ge $NMAX ]
then 
FINE="1"
fi
done
