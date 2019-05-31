for f in IR_*
do
YES=`tail -n 1 $f/veff*dat | awk '{if ($1<99999000000) print "1"; else print "0"}'`
if [ "$YES" == "1" ]
then
echo -n $f " "
fi
done
