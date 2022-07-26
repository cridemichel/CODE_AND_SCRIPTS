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
# $5 = jobfinished file (typically named '_DONE_')
# $6 ... possible extra args
import sys,os
args=sys.argv
jtype=int(args[1])
restart=args[2]
totsteps=int(args[3])
extsteps=int(args[4])
jobfinished=args[5]
### SIMULATION DEPENDENT VARIABLES
prepend='nohup mosrun ' # prepend to exec name
exec_name=prepend+'./mcsim_hardell' # sim executable
arg_start=' 2 ' # command line arg for starting
arg_restart=' 1 ' + restart #command line arg for restarting
##########################################################
# SIMULATION DEPENDENT CODE
# following functions depend on simulation file format
#
# get current simulation steps from restart file
def get_steps(restart):
    with open(restart) as f:
        firstline=f.readline()
        l=firstline.split(' ')
    return int(l[3])
#
#extend simulation by modifying restart file
def extend_sim(restart):   
    #os.system('rm '+ donefile)
    with open(restart) as f:
        ls=f.readlines()
    l=ls[1].split(' ')
    #overwrite file total steps with new value
    newsteps=int(l[1])+extsteps
    l[1]=str(newsteps)
    ls[1]=' '.join(l)
    with open(restart,"w") as f:
        for l in ls:
            f.write(l)
#
#establish whether simulation is finished (simulations steps > totsteps)
def sim_done(restart):
    if totsteps > 0:
        s=get_steps(restart)
        if s >= totsteps:
            return True
        else:
            return False
####################################################
if sim_done(restart):
    os.system('echo \"job regularly finished\" > '+jobfinished)
    quit()
if jtype==2:#start for the first time
    os.system(exec_name+arg_start + '> screen > 2>&1') 
else:#restart extending if requested, i.e. if totsteps > 0 
    if totsteps > 0:
        extend_sim(restart)
    os.system(exec_name + arg_restart + '>> screen 2>&1')
quit()
