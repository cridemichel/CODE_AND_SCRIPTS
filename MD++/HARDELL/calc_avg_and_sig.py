#!/usr/bin/env python3
import sys,os
import numpy as np
filen='_lista_dir_'
from operator import itemgetter
if (len(sys.argv)>1):
    ob=sys.argv[1]
else:
    ob='phi.dat'
os.system('ls -d P_* > ' + filen)
with open(filen) as f:
    lines=f.readlines()
lst=[]
lst2=[]
if len(sys.argv) > 2:
    fract=float(sys.argv[2].strip('\n'))
else:
    fract=0.9
for l in lines:
    os.chdir(l.strip('\n'))
    press = l.strip('\n').split('_')[-1]
    print('ob=', ob)
    with open(ob) as f:
        lob=f.readlines()
    #print('lphi=', lphi)
    nob=len(lob)
    snob=int(fract*nob)
    #print ('press=', press, 'nphi=', nphi, ' snphi=', snphi)
    #print('snphi=', snphi, ' fine=', nphi)
    if snob==nob:
        snob=nob-1
    subl=lob[snob:nob]
    #print('subl=', subl)
    sumob=0.0
    sumob2=0.0
    cc=0
    for ll in subl:
        obval=float(ll.strip('\n').split(' ')[1])
        #print('phi=',phi)
        sumob+=obval
        sumob2+=obval**2
        cc += 1.0
    lst.append([press,sumob/float(cc)])
    lst2.append([press,sumob2/float(cc)-(sumob/float(cc))**2])
    #print(sumphi/float(cc),' ',press)
    os.chdir('..')
lst.sort(key=itemgetter(0))
lst2.sort(key=itemgetter(0))
with open('avg.dat','w') as f:
    for l in lst:
        f.write(str(l[0])+ ' ' + str(l[1]) +'\n')
with open('stddev.dat','w') as f:
    for l in lst2:
        f.write(str(l[0])+ ' ' + str(l[1])+'\n')
os.system('rm ' + filen)
