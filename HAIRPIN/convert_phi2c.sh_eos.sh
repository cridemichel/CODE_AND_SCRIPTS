PHI2C="726.778"
FN="eos.dat"
if [ "$1" != "" ]
then
  FN="$1"
fi
cat $FN| awk -v k=$PHI2C '{print ($1*k,$2)}'
