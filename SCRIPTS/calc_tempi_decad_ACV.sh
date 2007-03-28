# 
echo "" > ./tempdecad_velACV_Phi$1.dat 
echo "" > ./tempdecad_omACV_Phi$1.dat
for f in X0_*
do
cd $f
if [ ! -d Phi$1 ]
then
cd ..
continue
fi
cd Phi$1
if [ ! -e velACV.dat ]
then 
cd ..
continue
fi
TEMPO=`../../find_first_below velACV.dat 0.367879`
echo "$f $TEMPO" >> ../../tempdecad_velACV_Phi$1.dat 
if [ ! -e omACV.dat ]
then 
cd ..
continue
fi
TEMPO=`../../find_first_below omACV.dat 0.367879` 
echo "$f $TEMPO" >> ../../tempdecad_omACV_Phi$1.dat
cd ..
cd ..
done
