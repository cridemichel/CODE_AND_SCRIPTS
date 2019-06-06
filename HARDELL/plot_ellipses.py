#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:24:44 2019

@author: demichel
"""
import math as m
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
#import matplotlib.transforms as transforms
import sys
nl=0
if len(sys.argv) > 1:
    filename=sys.argv[1]
else:
    filename="conf"
with open(filename) as f:
    for line in f:
        #print ("nl=", nl)
        if nl == 0:
            arr=line.strip('\n').split(' ')
            #print(arr[3])
            da=2.0*float(arr[3])
            db=2.0*float(arr[4])
            Lx=float(arr[1])
            Ly=float(arr[2])
            fig, ax = plt.subplots(figsize=(15, 15*Ly/Lx))
            plt.axis([-Lx/2.0-da, Lx/2.0+da, -Ly/2.0-db, Ly/2.0+db])
        else:
            arr=line.strip('\n').split(' ')
            x=(float(arr[0]),float(arr[1]))
            theta=180.0*np.arccos(float(arr[2]))/m.pi
            #print ("theta=", theta)
            ell=Ellipse(x,da,db,theta)
            ax.add_patch(ell)
        nl = nl + 1
plt.show()
