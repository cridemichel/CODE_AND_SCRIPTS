#!/usr/bin/env python3
import os,sys
os.system('ls -d X0_* > lista_dirs')
with open('lista_dirs') as f:
    lines=f.readlines()
for l in lines:
    os.chdir(l.strip('\n'))
    os.system('cp ../mcsim_hardell.sh .')
    os.system('../prepare_sim.py')
    os.system('cp res_X0_*.conf ~')
    os.chdir('..')
