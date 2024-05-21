for f in `ls -d X0_*`
do
	OD=`pwd`
	cd $f
	X0=$(echo $f | gawk -F _ '{print $2}')
	../prepare_all.py $X0
	cd $OD
done
