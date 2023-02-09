#!/usr/local/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import optimize
#import scipy.special as sf
doplot=False
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
def fgaussAD(th, a, norm):
    #ll=[]
    #K=np.sqrt(a)/(np.sqrt(pi)*sf.erfi(np.sqrt(a)))
    #for oneth in th:
    #    ll.append(K*np.exp(a*np.cos(oneth)*np.cos(oneth)))
    #return np.array(ll)
    #K=np.sqrt(a)/(np.sqrt(pi)*sf.erfi(np.sqrt(a)))
    return norm*np.exp(a*np.cos(th)*np.cos(th))
args=sys.argv
fi = 'distro.dat'
fiK = 'Kii.dat'
if len(args) > 2:
    rho=float(args[2])
if len(args) > 3:
    ii = int(args[3])
    if ii > 0:
        doplot=True
    else:
        doplot=False
if len(args) > 4:
    fi=args[4]
if len(args) > 5:
    fiK=args[5]
os.system('tail -n 1 ' + fiK +' > _aaa_')
with open('_aaa_',encoding='utf-8') as f:
    line=f.readlines()
l=line[0].strip('\n').split()
print(l)
usegauss=0
if len(args) > 1:
    usegauss=abs(int(args[1]))
os.system('rm _aaa_')
x_data, y_data = np.loadtxt(fi, usecols=(0,1), unpack=True)
if usegauss==0 or usegauss==3:
    print('Fitting Gauss function')
    params, params_covariance = optimize.curve_fit(fgauss, x_data, y_data, p0=[10.0])
    alpha=params[0]
    y_th = fgauss(x_data, alpha)
elif usegauss==1 or usegauss==4:
    print('Fitting Onsager function')
    params, params_covariance = optimize.curve_fit(fons, x_data, y_data, p0=[10.0])
    alpha=params[0]
    y_th = fons(x_data, alpha)
else:
    print('Fitting AD function')
    params, params_covariance = optimize.curve_fit(fgaussAD, x_data, y_data, p0=[1.0, 10.0])
    alpha=params[0]
    K=params[1]
    print('K=', K)
    y_th = fgaussAD(x_data, alpha, K)
print('alpha=', alpha, " params_covariance=", params_covariance)
if usegauss==0:
    fact=alpha * (1.0 - np.exp(-pi*pi*alpha/8.0))
    print('gauss fact=', fact)
elif usegauss==1:
    fact=alpha*np.tanh(alpha/2.0)
    print('onsager fact=', fact)
elif usegauss==2:
    fact=4.0*pi*(np.exp(alpha)-1.0)*K
    print('K=', K, ' AD fact=', fact)
else:
    print('Data from realigning simulation')
    fact=1.0
print('correcting factor with Parsons-Lee factor which is equal to', EtaP(veff*rho))
fact=fact*fact*EtaP(veff*rho)
K11=float(l[1])*fact
K22=float(l[2])*fact
K33=float(l[3])*fact
print(K11, ' ', K22, ' ', K33)
if doplot:
    plt.title(r'orientational distribution function')
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$f(\theta)$')
    plt.plot(x_data,y_data, 'xb',label='sim')
    plt.plot(x_data,y_th, '-y',label='fit')
    plt.savefig('angledistro.png')
    plt.legend()
    plt.show()
