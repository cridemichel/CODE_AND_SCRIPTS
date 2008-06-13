#PBS -S /bin/bash
#PBS -l walltime=1000000
#PBS -N _JOBNAME_

cp -r _LOCAL_FOLDER_ /scratch/fpiazza/$PBS_JOBID
cd /scratch/fpiazza/$PBS_JOBID
_EXE_NAME_ -fa ellipsoid_flex.par > screen 
