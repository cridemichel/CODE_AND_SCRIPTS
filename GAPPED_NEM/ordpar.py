#!/usr/bin/env python3
import sys,os
import numpy as np
filen='_lista_dir_'
from operator import itemgetter
os.system('ls -d P_*[^a-z]| sort -t _ -k 2 -n > ' + filen)
D=2.4
vol=np.pi*D*D/4.0*16.0
with open(filen) as f:
    lines=f.readlines()
lst=[]
if len(sys.argv) > 1:
    fract=float(sys.argv[1].strip('\n'))
else:
    fract=0.9
first=True
for l in lines:
    os.chdir(l.strip('\n'))
    print('dir=',l)
    press = l.strip('\n').split('_')[-1]
    with open('rho.dat') as f:
        lphi=f.readlines()
    #print('lphi=', lphi)
    nphi=len(lphi)
    snphi=int(fract*nphi)
    #print ('press=', press, 'nphi=', nphi, ' snphi=', snphi)
    #print('snphi=', snphi, ' fine=', nphi)
    if snphi==nphi:
        snphi=nphi-1
    subl=lphi[snphi:nphi]
    #print('subl=', subl)
    sumphi=0.0
    cc=0
    for ll in subl:
        phi=vol*float(ll.strip('\n').split(' ')[1])
        #print('phi=',phi)
        sumphi+=phi
        cc += 1.0
    os.system('ls cnf-*|sort -t - -k 2 -n | tail -n 10 > lista')
    os.system("../order_param lista | tail -n1 |awk '{print $3}' > _aaa_ ") 
    with open('_aaa_') as f:
        ll=f.readline().strip('\n')
    print('ll=',ll)
    lst.append([sumphi/float(cc),float(ll)])
    #print(sumphi/float(cc),' ',press)
    os.chdir('..')
#lst.sort(key=itemgetter(0))
with open('ordpar.dat','w') as f:
    for l in lst:
        f.write(str(l[0])+ ' ' + str(l[1]) +'\n')
os.system('rm ' + filen)
