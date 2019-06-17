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
if os.path.exists('lista_press'):
    with open('lista_press') as f:
        line=f.readline()
    lst2=[]
    lst=line.strip('\n').split(' ')        
    for l in lst:
        if l != '':
            lst2.append('P_'+l)
    with open('lista_dir','w') as f:
        for l in lst2:
            f.write(l+'\n')
else:    
    os.system('ls -d P_* > lista_dir')
with open("lista_dir") as f:
    lines=f.readlines()
for l in lines:
    os.chdir(l.strip('\n'))
    if os.path.exists('restart-0'):
        os.system('rm restart-[0,1]')
    if os.path.exists('cnf-final'):
        os.system('rm cnf-final')
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
    os.system('cp ../mcsim_hardell.sh .')
    os.chdir('..')
resfile='res_X0_' + X0 + '.conf'
os.system('echo \"-1 -1 60\" > '+resfile)
os.system('ls -d `pwd`/P_*/mcsim_hardell.sh >> '+resfile)
