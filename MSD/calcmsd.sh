PR=$HOME/ELLIPSOIDS/MSD/calcmsd
gunzip Store*gz
if [ "$1" == "" ]
then
NN=`cat Store-0-0.gz | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
else
NN=$1
fi
NN=`echo "$NN*5"| bc`
ls Store* | sort -t - -k 2 -k 3 -n > listamsd
$PR listamsd $NN 
gzip Store*
