#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.ticker import AutoMinorLocator
import matplotlib.transforms as mtransforms

matplotlib.rcParams['text.usetex'] = True
argv =sys.argv
#font = { 'size'   : 14}
#matplotlib.rc('font', **font)

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
#fig, axs = plt.subplots(2, 2, figsize=(16,9))
fig, axs = plt.subplot_mosaic([['(a)', '(b)'], ['(c)', '(d)']], figsize=(16,9), layout='constrained')

for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif')
#           bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))


#fig.tight_layout(pad=1.0)
#plt.subplots_adjust( wspace=0.1, hspace=0.1)

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
x, y = np.loadtxt('CASE_1_dt_0.001_4x7_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x2,y2= np.loadtxt('CASE_1_dt_0.001_4x7_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
sdf=2.0

ax = axs['(a)']
ax.plot(x, yp)
ax.plot(x, yp2)
#axs[0, 0].set_title(r'$dt = 0.001$') 
ax.set_ylim([0,4.5E-7])
ax.set(xlabel='', ylabel=r'$\Delta E$')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
#axs[0,0].plot(x, y, 'x',label='Traiettoria')
#
x, y = np.loadtxt('CASE_3_dt_0.003_4x7_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x2, y2 = np.loadtxt('CASE_3_dt_0.003_4x7_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
ax = axs['(b)']
ax.plot(x, yp)
ax.plot(x, yp2)
#axs[2, 0].set_ylim([0,3.0*1E-5*fact])
#axs[0, 1].set_title(r'$dt = 0.003$')
ax.set_ylim([0,6.0*1E-6])
ax.set(xlabel='', ylabel='')#ylabel=r'$\Delta E$')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.ticklabel_format(style='sci')
#axs[2, 1].plot(x, -y, 'tab:red')
x, y = np.loadtxt('CASE_5_dt_0.005_4x7_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x2, y2 = np.loadtxt('CASE_5_dt_0.005_4x7_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
ax=axs['(c)']
ax.plot(x, yp)
ax.plot(x, yp2)
#axs[2, 0].set_ylim([0,3.0*1E-5*fact])
#axs[1, 0].set_title(r'$dt = 0.005$')
ax.set_ylim([0,1.6*1E-5])
ax.set(xlabel='t (ns)', ylabel=r'$\Delta E$')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.ticklabel_format(style='sci')
#
x, y = np.loadtxt('CASE_10_dt_0.01_4x7_doublechain/delE.dat', comments=['#'], usecols=(0,1), unpack=True)
x2, y2 = np.loadtxt('CASE_10_dt_0.01_4x7_doublechain/runavg_delE.dat', comments=['#'], usecols=(0,1), unpack=True)
yp=np.array(y)
yp2=np.array(y2)
ax = axs['(d)']
ax.plot(x, yp)
ax.plot(x, yp2)
#axs[2, 0].set_ylim([0,3.0*1E-5*fact])
#axs[1, 1].set_title(r'$dt = 0.005$')
ax.set_ylim([0,8.5*1E-5])
ax.set(xlabel='t (ns)', ylabel='')#ylabel=r'$\Delta E$')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.ticklabel_format(style='sci')
plt.savefig('energy_case_SY.png')
plt.show()

