#!/usr/bin/env python3
import os,sys
arg=sys.argv
if len(arg)==1:
    os.system('ls -d P_* > lista_dirs')
with open('lista_dirs') as f:
    lines=f.readlines()
for l in lines:
    os.chdir(l.strip('\n'))
    os.system('cp ../job.py .')
    os.system('cp ../phcpars.asc .')
    os.system('cp ../mcsim_phc .')
    os.system('cp ../start.cnf .')
    os.system('../prepare_sim.py')
    os.chdir('..')
