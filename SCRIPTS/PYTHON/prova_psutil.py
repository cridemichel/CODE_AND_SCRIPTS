#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 18:09:50 2019

@author: demichel
"""
#import sys
import os
import psutil
pids = [pid for pid in psutil.pids()]
for pid in pids:
    uid=psutil.Process(pid).uids()[1]#use effective userid
    #print ('boh uid=',uid)
    #print ('uids=',psutil.Process(pid).uids())
   # print ('getuid=',os.getuid(),' uid=', uid)
    if int(uid) != int(os.getuid()):
        continue
    else:
        print ('cwd=', psutil.Process(pid).cwd())