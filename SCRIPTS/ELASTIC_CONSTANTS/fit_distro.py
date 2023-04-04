#!/usr/local/bin/python3
import matplotlib.pyplot as plt
import scipy.special as sf
import numpy as np
import sys
from scipy import optimize
doplot=False
pi=3.141592653589793
fi = 'distro.dat'
args=sys.argv
def fons(th, a):
    return (a/(4.0*pi*np.sinh(a))) * np.cosh(a*np.cos(th))
def fgauss(th, a):
    norm=np.sqrt(a)/(np.sqrt(pi)*sf.erfi(np.sqrt(a)))/(2.0*pi)
    return norm*np.exp(a*np.cos(th)*np.cos(th))
def fgaussct(cth, a):
    norm=np.sqrt(a)/(np.sqrt(pi)*sf.erfi(np.sqrt(a)))
    return norm*np.exp(a*cth*cth)
if len(args) > 2:
    ii = int(args[2])
    if ii > 0:
        doplot=True
    else:
        doplot=False
if len(args) > 3:
    fi=args[3]
usegauss=0
if len(args) > 1:
    usegauss=abs(int(args[1]))
x_data, y_data = np.loadtxt(fi, usecols=(0,1), unpack=True)
if usegauss==1:
    print('Fitting Onsager function')
    params, params_covariance = optimize.curve_fit(fons, x_data, y_data, p0=[10.0])
    alpha=params[0]
    y_th = fons(x_data, alpha)
elif usegauss==0:
    print('Fitting AD function')
    params, params_covariance = optimize.curve_fit(fgauss, x_data, y_data, p0=[1.0])
    alpha=params[0]
    y_th = fgauss(x_data, alpha)
else:
    print('Fitting AD (costh) function')
    params, params_covariance = optimize.curve_fit(fgaussct, x_data, y_data, p0=[1.0])
    alpha=params[0]
    y_th = fgaussct(x_data, alpha)
print('alpha=', alpha, " params_covariance=", params_covariance)
if doplot:
    plt.title(r'orientational distribution function')
    if usegauss > 1:
        plt.xlabel(r'$\cos(\theta)$')
        plt.ylabel(r'$f(\cos(\theta))$')
    else:
        plt.xlabel(r'$\theta$')
        plt.ylabel(r'$f(\theta)$')
    plt.plot(x_data,y_data, 'xb',label='sim')
    plt.plot(x_data,y_th, '-y',label='fit')
    plt.savefig('angledistro.png')
    plt.legend()
    plt.show()
