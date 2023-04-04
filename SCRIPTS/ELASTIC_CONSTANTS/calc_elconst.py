#!/usr/local/bin/python3
import os, sys
doplot=False
pi=3.141592653589793
veff=6.495
rho=0.0487
def EtaP(phi):
    return (4 - 3*phi)/(4*(1 - phi)*(1-phi))
args=sys.argv
fiK = 'Kii.dat'
if len(args) > 1:
    rho=float(args[1])
if len(args) > 2:
    fiK=args[2]
os.system('tail -n 1 ' + fiK +' > _aaa_')
with open('_aaa_',encoding='utf-8') as f:
    line=f.readlines()
l=line[0].strip('\n').split()
print(l)
usegauss=0
os.system('rm _aaa_')
print('correcting factor with Parsons-Lee factor which is equal to', EtaP(veff*rho))
fact=EtaP(veff*rho)
K11=float(l[1])*fact
K22=float(l[2])*fact
K33=float(l[3])*fact
print(K11, ' ', K22, ' ', K33)
