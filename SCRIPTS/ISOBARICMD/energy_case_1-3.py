#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.ticker import AutoMinorLocator

matplotlib.rcParams['text.usetex'] = True
argv =sys.argv
#font = { 'size'   : 14}
#matplotlib.rc('font', **font)

SMALL_SIZE = 22
MEDIUM_SIZE = 26
BIGGER_SIZE = 30
matplotlib.rc('font', size=SMALL_SIZE)          # controls default text sizes
matplotlib.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
matplotlib.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
matplotlib.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

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
fig, axs = plt.subplots(3, 2, figsize=(16,14))
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
xlbl=''#r'$t\;(\mathrm{ns})$'
ylbl='E'
logx=False
logy=False
yoshi='$\mathrm{SY+RESPA}$'
noyoshi='$\mathrm{NO\;\;SY+RESPA}$'
#plt.title(title)
plt.xlabel(xlbl)
t0=1E-3 # ns  
x, y = np.loadtxt('CASE_1_dt_0.001_4x7_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x=x*t0
x2,y2= np.loadtxt('CASE_1_dt_0.001_4x7_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
sdf=2.0
axs[0, 0].plot(x, yp, label=r'$E_{rel}$')
axs[0, 0].plot(x, yp2,label=r'$E_{cav}$')
axs[0, 0].set_title(yoshi)
#axs[0, 0].set_title(r'$\Delta t = 1\,\mathrm{fs}\;\;$',loc='right')
axs[0, 0].text(.95,.9,r'$\Delta t = 1\,\mathrm{fs}\;\;$', horizontalalignment='right', transform=axs[0,0].transAxes)
axs[0, 0].set_ylim([0,5.5E-7])
axs[0, 0].set(xlabel='', ylabel='') #ylabel=r'$\Delta E$')
axs[0, 0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0, 0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0, 0].legend()
#axs[0,0].plot(x, y, 'x',label='Traiettoria')
#
x, y = np.loadtxt('CASE_1_dt_0.001_noyoshi_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x=x*t0
x2,y2= np.loadtxt('CASE_1_dt_0.001_noyoshi_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
axs[0, 1].plot(x, yp)
axs[0, 1].plot(x, yp2)
axs[0, 1].set_title(noyoshi)
#axs[0, 1].set_title(r'$\Delta t = 1\,\mathrm{fs}\;\;$'+noyoshi)
axs[0, 1].set_ylim([0,5.5E-7])
axs[0, 1].xaxis.set_minor_locator(AutoMinorLocator())
axs[0, 1].yaxis.set_minor_locator(AutoMinorLocator())
axs[0, 1].set(xlabel='', ylabel='')
#axs[1, 0].plot(x, -y, 'tab:green')
###
x, y = np.loadtxt('CASE_2_dt_0.002_4x7_doublechain_3/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x=x*t0
x2,y2= np.loadtxt('CASE_2_dt_0.002_4x7_doublechain_3/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
axs[1, 0].plot(x, yp)
axs[1, 0].plot(x, yp2)
axs[1, 0].set_ylim([0,2.6*1E-6])
#axs[1, 0].set_title(r'$\Delta t = 2\,\mathrm{fs}\;\;$', loc='right')
axs[1, 0].text(.95,.9,r'$\Delta t = 2\,\mathrm{fs}\;\;$', horizontalalignment='right', transform=axs[1,0].transAxes)
axs[1, 0].xaxis.set_minor_locator(AutoMinorLocator())
axs[1, 0].yaxis.set_minor_locator(AutoMinorLocator())
axs[1, 0].set(xlabel='', ylabel='') #ylabel=r'$\Delta E$')
###
x, y = np.loadtxt('CASE_2_dt_0.002_noyoshi_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x=x*t0
x2, y2 = np.loadtxt('CASE_2_dt_0.002_noyoshi_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
axs[1, 1].plot(x, yp)
axs[1, 1].plot(x, yp2)
axs[1, 1].set_ylim([0,2.6*1E-6])
#axs[1, 1].set_title(r'$\Delta t = 2\,\mathrm{fs}\;\;$')
axs[1, 1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1, 1].yaxis.set_minor_locator(AutoMinorLocator())
axs[1, 1].set(xlabel='', ylabel='')
###
#axs[1, 1].plot(x, -y, 'tab:red')
#axs[2, 0].plot(x, -y, 'tab:green')
x, y = np.loadtxt('CASE_3_dt_0.003_4x7_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x=x*t0
x2, y2 = np.loadtxt('CASE_3_dt_0.003_4x7_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
axs[2, 0].plot(x, yp)
axs[2, 0].plot(x, yp2)
#axs[2, 0].set_ylim([0,3.0*1E-5*fact])
#axs[2, 0].set_title(r'$\Delta t = 3\,\mathrm{fs}\;\;$', loc='right')
axs[2, 0].text(.95,.9,r'$\Delta t = 3\,\mathrm{fs}\;\;$', horizontalalignment='right', transform=axs[2,0].transAxes)
axs[2, 0].set_ylim([0,7.2*1E-6])
axs[2, 0].set(xlabel=xlbl, ylabel='')#, ylabel=r'$\Delta E$')
axs[2, 0].xaxis.set_minor_locator(AutoMinorLocator())
axs[2, 0].yaxis.set_minor_locator(AutoMinorLocator())
axs[2, 0].ticklabel_format(style='sci')
#axs[2, 1].plot(x, -y, 'tab:red')
x, y = np.loadtxt('CASE_3_dt_0.003_noyoshi_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x=x*t0
x2,y2= np.loadtxt('CASE_3_dt_0.003_noyoshi_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
axs[2, 1].plot(x, yp)
axs[2, 1].plot(x, yp2)
axs[2, 1].set_ylim([0,7.2*1E-6])
#axs[2, 1].set_title(r'$\Delta t = 3\,\mathrm{fs}\;\;$'+noyoshi)
axs[2, 1].set(xlabel=xlbl, ylabel='')
axs[2, 1].xaxis.set_minor_locator(AutoMinorLocator())
axs[2, 1].yaxis.set_minor_locator(AutoMinorLocator())
plt.savefig('energy_case_1-3.png')
#for ax in axs.flat:
#    ax.set(xlabel='t', ylabel=r'$|(E(t)-E(0))/E(0)|$')
# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()
plt.show()

