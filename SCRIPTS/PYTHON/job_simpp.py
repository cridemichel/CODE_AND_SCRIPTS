#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 09:43:18 2019

@author: demichel
"""

#job script
# arguments passed by rescheduler.py
# $1 = 1 or 2; 1=restart 2=start
# $2 = restart file which is good for restarting 
# $3 = <totsteps> (-1 means to not extend)
# $4 = <extend steps> steps by which extend current
# $5 = jobfinished
# $5 ... possible extra args
import sys,os
args=sys.argv
type=int(args[1])
restart=args[2]
totsteps=int(args[3])
extsteps=int(args[4])
jobfinished=args[5]
exec_name='./mcsim_hardell'
def get_steps(restart):
    with open(restart) as f:
        firstline=f.readline()
        l=firstline.split(' ')
    #print ('l=',l)
    #print ('lines[0][2]=', ls[0][2])
    #overwrite filed total steps with new value
    return int(l[3])
def extend_sim(restart):   
    #os.system('rm '+ donefile)
    with open(restart) as f:
        ls=f.readlines()
    l=ls[0].split(' ')
    #print ('l=',l)
    #print ('lines[0][2]=', ls[0][2])
    #overwrite filed total steps with new value
    newsteps=int(l[2])+extsteps
    l[2]=str(newsteps)
    ls[0]=' '.join(l)
    #print ('dopo ls[0]=',ls[0])
    #print ('dir=', bn)
    #print ('res='+restart[w])
    with open(restart,"w") as f:
        for l in ls:
            f.write(l)
def sim_done(restart):
    if totsteps > 0:
        s=get_steps(restart)
        if s >= totsteps:
            return True
        else:
            return False
if sim_done(restart):
    os.system('echo job regularly finished > '+jobfinished)
    quit()
if type==2:#start for the first time
    os.system('nohup mosrun '+exec_name+' 1 > screen')
    quit()  
else:
    if totsteps > 0:
        extend_sim(restart)
    os.system('nohup mosrun '+exec_name + ' 2 ' + restart + ' > screen')
    quit()