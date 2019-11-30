#!/usr/bin/env python3
import os,sys
os.system('ls -d X0_* > lista_dirs')
with open('lista_dirs') as f:
    lines=f.readlines()
for l in lines:
    os.chdir(l.strip('\n'))
    os.system('ls -d P_* > _lista_dirs_')
    with open('_lista_dirs_') as f:
        linesP=f.readlines()
    for lP in linesP:
        d=lP.strip('\n')
        if os.path.exists(d):
            os.chdir(d)
            if os.path.exists('cnf-final'):
                os.system('cp cnf-final cnf-equilibrated')
            else:
                os.system('ls -d cnf-*[0-9] 2>&1 > _lista_files_')
                with open('_lista_files_') as f:
                    linesN=f.readlines()
                if len(linesN)==0:
                    print ('dir=', os.getcwd(), ' no cnf file found')
                    os.chdir('..')
                    continue
                ordlst=[]
                for lN in linesN:
                    lst=lN.strip('\n').split('-')
                    ordlst.append(lst[1])
                if len(ordlst) > 0:
                    ordlst.sort()
                    fn='cnf-'+ordlst[-1]
                    os.system('cp ' + fn + ' cnf-equilibrated')
            os.chdir('..')
    os.chdir('..')
