for d in `ls -d RUN_*`
do
OD=`pwd`
cd $d
FN=`ls dgeba_p_*`
nohup ./$FN -fa ellipsoid_dGEBA.par > screen & 
cd $OD
done
