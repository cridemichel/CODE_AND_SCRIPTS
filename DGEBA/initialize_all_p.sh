cd p_0.4701
../copy_files.sh Store-1-0
cd ..
cd p_0.5671
../copy_files.sh Store-2-0
cd ..
cd p_0.6313
../copy_files.sh Store-3-0
cd ..
cd p_0.6763
../copy_files.sh Store-4-0
cd ..
cd p_0.7096
../copy_files.sh Store-5-0
cd ..
cd p_0.7344
../copy_files.sh Store-6-0
cd ..
cd p_0.7682
../copy_files.sh Store-8-0
cd ..
cd p_0.8047
../copy_files.sh Store-12-0
cd ..
cd p_0.8248
../copy_files.sh Store-17-0
cd ..
for d in `ls -d p_*`
do 
  OD=`pwd`
  cd $d
  ../store_to_cnf_perm_bond.sh
  ../add_parameters_for_saving_presstens.sh
  for dd in `ls -d RUN_*`
  do
  ODD=`pwd`
  cd $dd
  ../../upd_param.sh storerate 5000000.0
  ../../upd_param.sh inifile start.cnf
  cp $HOME/UNIVSCAL_POLYMERS/autodgeba .
  ln -sf ./autodgeba dgeba_${d}_${dd} 
  cd $ODD
  done
  ../increase_Dt_all_runs.sh
  cd $OD
done
