PR=$HOME/ELLIPSOIDS/FQT/calcfqt
if [ "$1" == "" ]
then
NN=`cat Store-0-0.gz | gunzip -c | awk -F : '{if ($1=="NN") print $2}'` 
else
NN=$1
fi
echo "NN"=$NN
$PR 0  0  0  0 $NN 
$PR 2  0  2  0 $NN
$PR 0  0  1  0 $NN
$PR 0  0  2  0 $NN 
$PR 1 -1  1 -1 $NN 
$PR 1 -1  2 -1 $NN
$PR 1  0  1  0 $NN
$PR 1  0  2  0 $NN
$PR 1  1  1  1 $NN 
$PR 2  1  1  1 $NN
$PR 2 -2  2 -2 $NN
$PR 2  1  2  1 $NN
$PR 2  2  2  2 $NN
