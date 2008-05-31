for fff in NVE-T*
do
cd $fff
LLF=`tail -1 CnfT*-NVE-1`
for ffff in Cnf*
do
LL=`tail -1 $ffff`
if [ "$LL" != "$LLF" ]
then
echo $fff/$ffff "is a corrupted file..." 
rm -f $ffff
fi
done
cd ..
done
