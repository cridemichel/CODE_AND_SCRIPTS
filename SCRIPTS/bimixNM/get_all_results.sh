DATA=`date '+%d%m%y'`
FL=""
for f in NM*-*
do
FL="$FL $f/rho*/NVE-T*/Fqs* $f/rho*/NVE-T*/N-sqt* $f/rho*/NVE-T*/MSD* $f/rho*/NVE-T*/Sq* $f/rho*/NVE-T*/alpha*" 
FL="$FL $f/rho*/T*/Fqs* $f/rho*/T*/N-sqt* $f/rho*/T*/MSD* $f/rho*/T*/Sq* $f/rho*/T*/alpha*"
done
echo "FILE_LIST=" $FL
tar czf  BMNMres${DATA}.tgz $FL
#tar czf BMNMres${DATA}-2.tgz NM11-6/rho*/NVE-T*/Fqs* NM11-6/rho*/NVE-T*/N-sqt* NM11-6/rho*/NVE-T*/MSD* NM11-6/rho*/NVE-T*/Sq* NM11-6/rho*/NVE-T*/alpha* NM12-6/rho*/NVE-T*/Fqs* NM12-6/rho*/NVE-T*/N-sqt* NM12-6/rho*/NVE-T*/MSD* NM12-6/rho*/NVE-T*/Sq* NM12-6/rho*/NVE-T*/alpha* NM9-6/rho*/NVE-T*/Fqs* NM9-6/rho*/NVE-T*/N-sqt* NM9-6/rho*/NVE-T*/MSD* NM9-6/rho*/NVE-T*/Sq* NM9-6/rho*/NVE-T*/alpha* NM12-11/rho*/NVE-T*/Fqs* NM12-11/rho*/NVE-T*/N-sqt* NM12-11/rho*/NVE-T*/MSD* NM12-11/rho*/NVE-T*/Sq* NM12-11/rho*/NVE-T*/alpha* NM8-5/rho*/NVE-T*/Fqs* NM8-5/rho*/NVE-T*/N-sqt* NM8-5/rho*/NVE-T*/MSD* NM8-5/rho*/NVE-T*/Sq* NM8-5/rho*/NVE-T*/alpha*
