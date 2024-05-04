#!/bin/bash
for s in `ls submit_p_*`
do
sbatch ./$s
done
