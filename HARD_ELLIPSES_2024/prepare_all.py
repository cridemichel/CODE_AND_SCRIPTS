#!/usr/bin/env python3
import os,sys
arg=sys.argv
if len(arg)==2:
    os.system('ls -d P_* > lista_dirs')
with open('lista_dirs') as f:
    lines=f.readlines()
for l in lines:
    XDIR=sys.argv[1]
    os.chdir(l.strip('\n'))
    PDIR=l.strip('\n')
    os.system('cp ../../job.py .')
    os.system('cp ../../hellpars_templ.xasc hellpars.xasc')
    os.system('cp ../../mcsim_hardellipses .')
    os.system('cp ../start.cnf .')
    os.system('sh ../../upd_parfile.sh ' + XDIR + ' ' + PDIR)
    os.chdir('..')
