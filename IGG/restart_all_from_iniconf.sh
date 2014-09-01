for f in `cat $1`
do
cd $f 
EXE=`ls ellips-IgG*`
EXEN="N-$EXE"
rm COORD_TMP[0,1]
cat COORD_TMP_ASCII0| awk 'BEGIN{skip=0} {if (doout==1 && $1!="RF" && skip==0 && $1!="saveBonds:" && $1!="maxbondsSaved:" && $1!="nintersIJ:") print $0; if ($0=="1 1 4 0 1 0 0 1") skip=1; if ($1=="@@@") {at++; if (skip==1) {skip=0; print ("@@@");}}; if ($1=="@@@" && at==1) doout=1;}' > CorIni
ln -sf ./ellipsoid $EXEN 
echo "EXEN=" $EXEN
nohup mosrun ./$EXEN -fa ellipsoid_flex.par > screen &
sleep 0.5
cd ..
done
