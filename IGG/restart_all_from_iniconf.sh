for f in `cat $1`
do
cd $f 
EXE=`ls ellips-IgG*`
EXEN="N-$EXE"
cp $HOME/IGG/EDBD/ellipsoid-new .
rm COORD_TMP[0,1] 2> /dev/null
rm Store-*-* 2> /dev/null
WC0=`wc -l COORD_TMP_ASCII0|awk '{print $1}'`
WC1=`wc -l COORD_TMP_ASCII1|awk '{print $1}'`
if [  $WC0 -gt $WC1 ]
then
COF="COORD_TMP_ASCII0"
else
COF="COORD_TMP_ASCII1"
fi
cat $COF| awk 'BEGIN{skip=0} {if (doout==1 && $1!="RF" && skip==0 && $1!="saveBonds:" && $1!="maxbondsSaved:" && $1!="nintersIJ:") print $0; if ($0=="1 1 4 0 1 0 0 1") skip=1; if ($1=="@@@") {at++; if (skip==1) {skip=0; print ("@@@");}}; if ($1=="@@@" && at==1) doout=1;}' > CorIni
ln -sf ./ellipsoid-new $EXEN 
echo "EXEN=" $EXEN
../../set_params.sh ellipsoid_flex.par seed -1 intervalSum 5000.0
nohup mosrun ./$EXEN -fa ellipsoid_flex.par > screen &
sleep 0.5
cd ..
done
