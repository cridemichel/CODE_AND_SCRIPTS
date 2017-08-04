# 2 = K22; 3 = K11; 4 = k33  
# $1 = ALPHA
if [ "$1" == "" ]
then
echo "calc_vexcl,sh <alpha> <type=1,2,3 1=K11, 2=K22, 3=K33>"
exit
fi
if [ "$1" != "" ]
then
ALPHA="$1"
else
ALPHA="20.0"
fi
if [ "$2" != "" ]
then
if [ "$2" == "1" ]
then
TYPE="3"
elif [ "$2" == "2" ]
then
TYPE="2"
elif [ "$2" == "3" ]
then	
TYPE="4"
fi
else 
TYPE="3"
fi
echo "ALPHA= " $ALPHA "TYPE= " $TYPE
./CG-CHROMONICS-ELCONST-MC 1.1 10 $ALPHA 20000000000 $TYPE 10000 
