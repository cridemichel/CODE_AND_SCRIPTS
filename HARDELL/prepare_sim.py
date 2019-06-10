#!/usr/bin/env python3
import os, sys
os.system('../create_dirs.py')
press='1.0'
curd=os.getcwd()
d=os.path.split(curd.strip('\n'))[1]
X0=d.split('_')[1].strip(' ')
sa=str(float(X0)*0.5)
print ('X0=',X0)
os.system('cp ../hepars_templ.asc .')
lines2new=[]
with open("hepars_templ.asc") as f:
    lines2=f.readlines()
    for l2 in lines2:
        #print ('l2=', l2)
        ls=l2.strip('\n').split(' ')
        ln = ' '.join([l.replace('_PRESS_',str(press)) for l in ls])+'\n'
        ln = ln.replace('_A_',sa)
        lines2new.append(ln)
with open("hepars.asc","w") as f:
    for l in lines2new:
        #print ('l=',l)
        f.write(l)
os.system('../create_conf.sh')
os.system('ls -d P_* > lista_dir')
with open("lista_dir") as f:
    lines=f.readlines()
for l in lines:
    os.chdir(l.strip('\n'))
    ll=l.strip('\n').split('_')
    press=ll[1]
    lines2new=[]
    os.system('cp ../hepars_templ.asc .')
    with open("hepars_templ.asc") as f:
        lines2=f.readlines()
        for l2 in lines2:
            #print ('l2=', l2)
            ls=l2.strip('\n').split(' ')
            ln = ' '.join([l.replace('_PRESS_',str(press)) for l in ls])+'\n'
            ln = ln.replace('_A_',sa)
            lines2new.append(ln)
    with open("hepars.asc","w") as f:
        for l in lines2new:
            #print ('l=',l)
            f.write(l)
    os.system('cp ../startcnf .')
    os.system('cp ../mcsim_hardell .')
    os.chdir('..')
resfile='res_X0_' + X0 + '.conf'
os.system('echo \"-1 -1 10\" > '+resfile)
os.system('ls -d `pwd`/P_*/*mcsim* >> '+resfile)
