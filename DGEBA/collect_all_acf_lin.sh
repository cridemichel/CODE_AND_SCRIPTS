if [ "$1" = "" ]
then
DT=1.0
else
DT=$1
fi
cd p_0.8248
../calc_all_acf.sh
../media_dbl_acf_lin.exe lista_acf $DT > media_acf_lin.dat
cd ..
cd p_0.8047
../calc_all_acf.sh
../media_dbl_acf_lin.exe lista_acf $DT > media_acf_lin.dat
cd ..
cd p_0.7682
../calc_all_acf.sh
../media_dbl_acf_lin.exe lista_acf $DT > media_acf_lin.dat
cd ..
