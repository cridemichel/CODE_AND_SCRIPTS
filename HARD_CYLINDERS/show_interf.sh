if [ "$2" != "" ]
then
THRA=$2
else
THRA="0.0"
fi
if [ "$3" != "" ]
then
THRB=$3
else
THRB="44.0"
fi
LZ="94.0"
COLA="turquoise"
COLB="green"
cat $1 | gawk -v lz="$LZ" -v cola=$COLA -v colb=$COLB -v thra="$THRA" -v thrb="$THRB" '{if (NF == 10 && NR!=1) {if (($3 < thrb) && ($3 > thra)) {for (i=1;i<=9;i++) printf("%s ",$i); printf(" C[%s]\n",cola);} else {for (i=1;i<=9;i++) printf("%s ",$i); printf("C[%s]\n",colb);}}; if (NF < 10) {print $0;};}'
