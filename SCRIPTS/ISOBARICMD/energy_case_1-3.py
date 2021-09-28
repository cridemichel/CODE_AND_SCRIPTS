#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

matplotlib.rcParams['text.usetex'] = True
argv =sys.argv
absenedr=True
def func(xx, aa, bb):
    return aa * xx  + bb

def calcstddev(vec):
    m=0
    c1=0
    for v in vec:
        m += v
        c1+=1
    m = m/c1
    s=0
    for v in vec:
        s += (v-m)*(v-m)
    s = np.sqrt(s / c1)
    return s

def calcedr(vec,lst,lst2):
    # calc energy drift according to
    #Yu et al. Chemical Physics 370 (2010) 294–305 (see Fig. 1 therein)
    c1=1.0
    v0 = vec[0]
    ssum=0.0
    for v in vec:
        if absenedr==True:
            val=np.abs((v-v0)/v0)
        else:
            val=(v-v0)/v0
        ssum += val
        lst.append(ssum/c1)
        lst2.append(val)
        c1+=1.0

def calcra(vec, lst):
    c1=1.0
    ssum=0.0
    for v in vec:
        ssum += v
        lst.append(ssum/c1)
        c1+=1.0

#if len(argv) < 2:
#    print("Please supply a file with data to plot!")
#    quit()
#parse arguments
ra=[1.0]
ra.clear()
runavg=False
dofit=False
calcenedrift=False
fig, axs = plt.subplots(3, 2)
#  def plot_addons():
#      calcedr(y, ra, yp)
#      sd=calcstddev(y)
#      if dofit:
#              popt, pcov = curve_fit(func, x, y)
#      av=np.mean(y)
#      sd=np.std(y)
#      plt.ylim(av-sd*sdf,av+sd*sdf)
#      if runavg:
#          calcra(y,yp2)
#          plt.plot(x, yp2,'-g', linewidth=2.0)
#      if calcenedrift:
#         plt.plot(x, ra,'-r', linewidth=2.0)
#      if dofit:
#         plt.plot(x, func(x, *popt), '-', label='fit: y=(%G)*x + (%G)' % tuple(popt))
#y=yp
#ax.tick_params(labeltop=True, labelright=True)
title=''
xlbl='t'
ylbl='E'
logx=False
logy=False
#plt.title(title)
plt.xlabel(xlbl)
x, y = np.loadtxt('CASE_1_dt_0.001_4x7_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[1.0]
yp.clear()
calcedr(y, ra, yp)
fact=1E7
y=fact*np.array(yp)
sd=calcstddev(y)
yp2=[1.0]
yp2.clear()
calcra(y,yp2)
sdf=2.0
axs[0, 0].plot(x, y)
axs[0, 0].plot(x, yp2)
axs[0, 0].text(0, 6E-7*fact, r'$dt = 0.001\; n_{ys}=4\; n_r = 7$')
axs[0, 0].set_ylim([0,8E-7*fact])
axs[0, 0].set(xlabel='', ylabel=r'$\Delta E\; \times 10^{-7}$')
#axs[0,0].plot(x, y, 'x',label='Traiettoria')
#
x, y = np.loadtxt('CASE_1_dt_0.001_noyoshi_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[]
calcedr(y, ra, yp)
fact=1E7
y=fact*np.array(yp)
sd=calcstddev(y)
yp2=[]
calcra(y,yp2)
axs[0, 1].plot(x, y)
axs[0, 1].plot(x, yp2)
axs[0, 1].text(0, 6E-7*fact, r'$dt = 0.001\; n_{ys}=1\; n_r = 1$')
axs[0, 1].set_ylim([0,8E-7*fact])
axs[0, 1].set(xlabel='', ylabel='')
#axs[1, 0].plot(x, -y, 'tab:green')
###
x, y = np.loadtxt('CASE_2_dt_0.002_4x7_doublechain_3/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[]
calcedr(y, ra, yp)
fact=1E6
y=fact*np.array(yp)
sd=calcstddev(y)
yp2=[]
calcra(y,yp2)
axs[1, 0].plot(x, y)
axs[1, 0].plot(x, yp2)
axs[1, 0].set_ylim([0,3.0*1E-6*fact])
axs[1, 0].set_title(r'$dt = 0.002\; n_{ys}=4\; n_r = 7$')
axs[1, 0].set(xlabel='', ylabel=r'$\Delta E \times 10^{-6}$')
###
x, y = np.loadtxt('CASE_2_dt_0.002_noyoshi_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[]
calcedr(y, ra, yp)
y=fact*np.array(yp)
sd=calcstddev(y)
yp2=[]
calcra(y,yp2)
axs[1, 1].plot(x, y)
axs[1, 1].plot(x, yp2)
axs[1, 1].set_ylim([0,3.0*1E-6*fact])
axs[1, 1].set_title(r'$dt = 0.002\; n_{ys}=1\; n_r = 1$')
axs[1, 1].set(xlabel='', ylabel='')
###
#axs[1, 1].plot(x, -y, 'tab:red')
#axs[2, 0].plot(x, -y, 'tab:green')
x, y = np.loadtxt('CASE_3_dt_0.003_4x7_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[]
calcedr(y, ra, yp)
y=yp
sd=calcstddev(y)
yp2=[]
calcra(y,yp2)
axs[2, 0].plot(x, y)
axs[2, 0].plot(x, yp2)
axs[2, 0].set_title(r'$dt = 0.003\; n_{ys}=4\; n_r = 7$')
axs[2, 0].set(xlabel='t (ns)', ylabel=r'$|(E(t)-E(0))/E(0)|$')
#axs[2, 1].plot(x, -y, 'tab:red')
x, y = np.loadtxt('CASE_3_dt_0.003_noyoshi_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[]
calcedr(y, ra, yp)
y=yp
sd=calcstddev(y)
yp2=[]
calcra(y,yp2)
axs[2, 1].plot(x, y)
axs[2, 1].plot(x, yp2)
axs[2, 1].set_title(r'$dt = 0.003\; n_{ys}=1\; n_r = 1$')
axs[2, 1].set(xlabel='t (ns)', ylabel='')
plt.savefig('energy_case_1-3.png')
#for ax in axs.flat:
#    ax.set(xlabel='t', ylabel=r'$|(E(t)-E(0))/E(0)|$')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()
plt.show()

