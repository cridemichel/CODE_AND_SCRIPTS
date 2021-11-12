import numpy as np
import sys
#Usage calc_folding.py <threshold_in_degree> 
with open('cnf-final',encoding='utf-8') as f:
    lines=f.readlines()
thr=float(sys.argv[1])
lines=lines[1:]
orients=[]
for l in lines:
    ls=l.strip('\n').split()
    orients.append(np.array([float(ls[3]),float(ls[4]),float(ls[5])]))
Np=int(len(orients)/2)
nfolded=0
for i in range(0,Np):
    i1=i*2
    i2=i*2+1
    theta=180.0*np.arccos(np.dot(orients[i1],orients[i2]))/np.pi
    if theta < thr: # if theta is less than 10 degree
        nfolded+=1
print('# of gapped duplexes foldes is ', nfolded)
