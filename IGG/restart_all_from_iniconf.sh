for f in `cat $1`
do
cd $f 
EXE=`ls ellips-IgG*`
cat COORD_TMP_ASCII0| awk 'BEGIN{skip=0} {if (doout==1 && $1!="RF" && skip==0) print $0;  if ($0=="1 1 4 0 1 0 0 1") skip=1; if ($1=="@@@") {at++; if (skip==1) {skip=0; print ("@@@");}}; if ($1=="@@@" && at==1) doout=1;}' > CorIni
nohup mosrun ./$EXE -fa ellipsoid_flex.par > screen &
sleep 0.5
cd ..
done
