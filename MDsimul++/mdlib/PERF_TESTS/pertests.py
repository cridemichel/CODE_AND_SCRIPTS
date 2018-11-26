#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 17:04:26 2018

@author: demichel
"""
import matplotlib as mp
import numpy as np  
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pylab as plt
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
mp.rc('font', **font)  
minorLocatœorx = MultipleLocator(5)
majorLocatorx = MultipleLocator(50)
minorLocatory = MultipleLocator(0.2)
majorLocatory = MultipleLocator(1)
su=np.loadtxt("speedup_ALL_3.out","double")
plt.title(r'$Test$')
plt.xlabel("matrix order",fontsize=15)
plt.ylabel("speedup", fontsize=15)
#fig, ax = plt.subplots()
#fig, ax = plt.subplots() 
#ax.xaxis.set_major_locator(majorLocatorx)
#ax.xaxis.set_minor_locator(minorLocatorx)
#ax.yaxis.set_minor_locator(minorLocatory)
#ax.yaxis.set_major_locator(majorLocatory)
arma, =plt.plot(su[:,[0]],su[:,[1]],'r--',label='Armadillo',marker=8) # : = all rows and [1] means first column
eigen,=plt.plot(su[:,[0]],su[:,[2]],'b--',label='Eigen',marker='^')
plt.legend(handles=[arma, eigen])
plt.axis([0, 200, 0, 8])
plt.savefig('miscops_test.pdf',format='pdf',dpi=600)