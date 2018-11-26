#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 17:04:26 2018

@author: demichel
"""
import matplotlib as mp
import numpy as np  
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.pylab as plt
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
mp.rc('font', **font)
minorLocatorx = MultipleLocator(5) # define axes major and minor ticks
majorLocatorx = MultipleLocator(50)
minorLocatory = MultipleLocator(0.2)
majorLocatory = MultipleLocator(1)
su=np.loadtxt("speedup_ALL_3.out","double")
gl=np.loadtxt("speedup_glm.out", "double")
fig, ax = plt.subplots()  # get axes of plots
plt.title(r'$Test\;\; M_1 = M_1 + M_1 M_2 M_4^{-1} + (M_1 + M_2 + M_3 + M_4)$',fontsize=12) # set title enabling latex parsing
plt.xlabel("matrix order",fontsize=15) # xaxis label with font size
plt.ylabel("speedup", fontsize=15)      # yaxis label with font size
plt.xticks(np.arange(0, 200, step=20))  # sey major ticks with labels
plt.yticks(np.arange(0, 10, step=1))   #   set major ticks with labels
ax.xaxis.set_minor_locator(minorLocatorx)#set minor ticks of xaxis
ax.yaxis.set_minor_locator(minorLocatory)#set minor ticks of yaxis
arma, =plt.plot(su[:,[0]],su[:,[1]],'r--',label='Armadillo',marker=8) # : = all rows and [1] means first column
eigen,=plt.plot(su[:,[0]],su[:,[2]],'b--',label='Eigen',marker='^')
plt.legend(handles=[arma, eigen])
plt.axis([0, 200, 0, 8])
plt.savefig('miscops_test.pdf',format='pdf',dpi=600)
plt.show()