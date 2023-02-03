import numpy as np
import os
from scipy import optimize
pi=3.141592653589793
def fons(th, a):
    return (a/(4.0*pi*np.sinh(a))) * np.cosh(a*np.cos(th))
os.system('tail -n 1 Kii.dat > _aaa_')
with open('_aaa_') as f:
    line=f.readlines()
l=line[0].strip('\n').split()
print(l)
x_data, y_data = np.loadtxt('distro.dat', usecols=(0,1), unpack=True)
params, params_covariance = optimize.curve_fit(fons, x_data, y_data, p0=[10.0])
alpha=params[0]
fact=alpha*np.tanh(alpha/2.0)
fact=fact*fact
K11=float(l[1])*fact
K22=float(l[2])*fact
K33=float(l[3])*fact
print(K11, ' ', K22, ' ', K33)

