#PBS -S /bin/bash
#PBS -l walltime=1000000
#PBS -N _JOBNAME_

cd _LOCAL_FOLDER_
echo _JOBNAME_ $PBS_JOBID > JOB_ID.txt
./ellipsoid -c >> screen 
