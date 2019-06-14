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
# $4 = <extend steps> steps by which extend current job
# $5 = jobfinished file (tyipically '_DONE_')
# $5 ... possible extra args
import sys,os
args=sys.argv
type=int(args[1])
restart=args[2]
totsteps=int(args[3])
extsteps=int(args[4])
jobfinished=args[5]
prepend='nohup mosrun '
exec_name=prepend+'./mcsim_hardell'
arg_start=' 2 '
arg_restart=' 1 ' + restart
def get_steps(restart):
    with open(restart) as f:
        firstline=f.readline()
        l=firstline.split(' ')
    return int(l[3])
def extend_sim(restart):   
    #os.system('rm '+ donefile)
    with open(restart) as f:
        ls=f.readlines()
    l=ls[0].split(' ')
    #overwrite filed total steps with new value
    newsteps=int(l[2])+extsteps
    l[2]=str(newsteps)
    ls[0]=' '.join(l)
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
    os.system('echo \"job regularly finished\" > '+jobfinished)
    quit()
if type==2:#start for the first time
    os.system(exec_name+arg_start + ' > screen') 
else:#restart extending if requested, i.e. if totsteps > 0 
    if totsteps > 0:
        extend_sim(restart)
    os.system(exec_name + arg_restart + ' > screen')
quit()