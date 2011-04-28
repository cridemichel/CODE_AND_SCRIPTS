MODEL="2"
TEMP="0.12"
SIGMA="1.21667" # sporge di 0.15 ed ha quindi un diametro esterno di 0.4 D come richiesto da Tommaso
PHIE="0.02"
PHI="0.02"
FINISHED="0"
CC="1"
while [ "$FINISHED" == "0" ]
do 
EQSTPS=`echo "2000 $CC"| awk '{printf("%d",$1*$2)}'`
./sim1statepnt_COKECAN_MCNVT.sh $PHI 10000000 2.0 $TEMP $SIGMA $MODEL 50000 &
sleep 2
PHI=`echo "$PHI 0.02"| awk '{printf("%1.2f",$1+$2)}'`
FINISHED=`echo $PHI | awk -v phi=$PHI -v phie=$PHIE '{if (phi > phie) printf("1"); else printf("0");}'`
CC=$[$CC+1]
done
