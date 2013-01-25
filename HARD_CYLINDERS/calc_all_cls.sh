ls -d N_* > listadir_
for f in `cat listadir_`
do
cd $f
ls Cnf*  > listafiles
../../newnclust-due-many-sq
ls NdcCnf* > listacls
../../calc_nu_avg listacls > Ndc_avg.dat
cd ..
done
rm listadir_
