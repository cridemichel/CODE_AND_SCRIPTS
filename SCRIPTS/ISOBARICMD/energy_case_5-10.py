#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.ticker import AutoMinorLocator

matplotlib.rcParams['text.usetex'] = True
argv =sys.argv
absenedr=True
def func(xx, aa, bb):
    return aa * xx  + bb

SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 24

matplotlib.rc('font', size=SMALL_SIZE)          # controls default text sizes
matplotlib.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
matplotlib.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
matplotlib.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

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
fig, axs = plt.subplots(2, 2, figsize=(16,9.4))
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
plt.title(title)
plt.xlabel(xlbl)
x, y = np.loadtxt('CASE_5_dt_0.005_4x7_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[1.0]
yp.clear()
calcedr(y, ra, yp)
y=yp
sd=calcstddev(y)
yp2=[1.0]
yp2.clear()
calcra(y,yp2)
sdf=2.0
axs[0, 0].plot(x, y)
axs[0, 0].plot(x, yp2)
axs[0, 0].set_title(r'$dt = 0.005\; n_{ys}=4\; n_r = 7$')
axs[0, 0].set_ylim([0,3E-5])
axs[0, 0].set(xlabel='', ylabel=r'$\Delta E$')
axs[0, 0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0, 0].yaxis.set_minor_locator(AutoMinorLocator())
#axs[0,0].plot(x, y, 'x',label='Traiettoria')
#
x, y = np.loadtxt('CASE_5_dt_0.005_noyoshi_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[]
calcedr(y, ra, yp)
y=yp
sd=calcstddev(y)
yp2=[]
calcra(y,yp2)
axs[0, 1].plot(x, y)
axs[0, 1].plot(x, yp2)
axs[0, 1].set_title(r'$dt = 0.005\; n_{ys}=1\; n_r = 1$')
axs[0, 1].set_ylim([0,3E-5])
axs[0, 1].set(xlabel='', ylabel='')
axs[0, 1].xaxis.set_minor_locator(AutoMinorLocator())
axs[0, 1].yaxis.set_minor_locator(AutoMinorLocator())
#axs[1, 0].plot(x, -y, 'tab:green')
###
x, y = np.loadtxt('CASE_10_dt_0.01_4x7_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[]
calcedr(y, ra, yp)
y=yp
sd=calcstddev(y)
yp2=[]
calcra(y,yp2)
axs[1, 0].plot(x, y)
axs[1, 0].plot(x, yp2)
axs[1, 0].set_ylim([0,1.25E-4])
axs[1, 0].set_title(r'$dt = 0.01\; n_{ys}=4\; n_r = 7$')
axs[1, 0].set(xlabel=r'$t (ns)$', ylabel=r'$\Delta E$')
axs[1, 0].ticklabel_format(axis='y', style='sci',scilimits=(0,0))
axs[1, 0].xaxis.set_minor_locator(AutoMinorLocator())
axs[1, 0].yaxis.set_minor_locator(AutoMinorLocator())

###
x, y = np.loadtxt('CASE_10_dt_0.01_noyoshi_doublechain/energy.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=[]
calcedr(y, ra, yp)
y=yp
sd=calcstddev(y)
yp2=[]
calcra(y,yp2)
axs[1, 1].plot(x, y)
axs[1, 1].plot(x, yp2)
axs[1, 1].set_ylim([0,1.25E-4])
axs[1, 1].set_title(r'$dt = 0.01\; n_{ys}=1\; n_r = 1$')
axs[1, 1].ticklabel_format(axis='y', style='sci',scilimits=(0,0))
axs[1, 1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1, 1].yaxis.set_minor_locator(AutoMinorLocator())

#axs[1, 1].set(xlabel=r'$t (ns)$', ylabel='')
###
#axs[1, 1].plot(x, -y, 'tab:red')
axs[1, 1].set_title(r'$dt = 0.01\; n_{ys}=1\; n_r = 1$')
#axs[2, 0].plot(x, -y, 'tab:green')
plt.savefig('energy_case_5-10.png')
#for ax in axs.flat:
#    ax.set(xlabel='t', ylabel=r'$|(E(t)-E(0))/E(0)|$')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()
plt.show()

