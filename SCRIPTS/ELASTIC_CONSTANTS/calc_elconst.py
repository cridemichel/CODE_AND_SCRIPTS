#!/usr/local/bin/python3
import numpy as np
import os, sys
from scipy import optimize
pi=3.141592653589793
veff=6.495
rho=0.0487
def fons(th, a):
    return (a/(4.0*pi*np.sinh(a))) * np.cosh(a*np.cos(th))
def EtaP(phi):
    return (4 - 3*phi)/(4*(1 - phi)*(1-phi))
def fgauss(th, a):
    ll=[]
    for oneth in th:
        if oneth <= pi*0.5:
            arg = oneth*oneth
        else:
            arg = (pi-oneth)*(pi-oneth)
        ll.append((a/(4*pi))*np.exp(-a*arg*0.5))
    return np.array(ll)
args=sys.argv
fi = 'distro.dat'
fiK = 'Kii.dat'
if len(args) > 2:
    rho=float(args[2])
if len(args) > 3:
    fi=args[3]
if len(args) > 4:
    fiK=args[4]
os.system('tail -n 1 ' + fiK +' > _aaa_')
with open('_aaa_',encoding='utf-8') as f:
    line=f.readlines()
l=line[0].strip('\n').split()
print(l)
usegauss=0
if len(args) > 1:
    usegauss=int(args[1])
os.system('rm _aaa_')
x_data, y_data = np.loadtxt(fi, usecols=(0,1), unpack=True)
if usegauss==0:
    params, params_covariance = optimize.curve_fit(fgauss, x_data, y_data, p0=[10.0])
    alpha=params[0]
    print('alpha=', alpha)
elif usegauss==1:
    params, params_covariance = optimize.curve_fit(fons, x_data, y_data, p0=[10.0])
    alpha=params[0]
    print('alpha=', alpha)
if usegauss==0:
    fact=alpha * (1.0 - np.exp(-pi*pi*alpha/8.0))
elif usegauss==1:
    fact=alpha*np.tanh(alpha/2.0)
else:
    fact=1.0
fact=fact*fact*EtaP(veff*rho)
K11=float(l[1])*fact
K22=float(l[2])*fact
K33=float(l[3])*fact
print(K11, ' ', K22, ' ', K33)
