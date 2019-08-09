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
# $5 ... possible extra args
import sys,os
args=sys.argv
jtype=int(args[1])
restart=args[2]
totsteps=int(args[3])
extsteps=int(args[4])
jobfinished=args[5]
### SIMULATION DEPENDENT VARIABLES
prepend='nohup mosrun ' # prepend to exec name
exec_name=prepend+'./ellipsoid ' # sim executable
arg_start=' -fa ellipsoid_flex.par ' # command line arg for starting
arg_restart=' -ca ' + restart #command line arg for restarting
##########################################################
# SIMULATION DEPENDENT CODE
# following functions depend on simulation file format
#
# get current simulation steps from restart file
def get_steps(restart):
    with open(restart) as f:
        lines=f.readlines
    for l in lines:
        par=l.strip('\n').split(':').strip(' ')
        if par[0] == 'curStep':
            stps=int(par[1])
            break
    return stps
#
#extend simulation by modifying restart file
def extend_sim(restart):   
    #os.system('rm '+ donefile)
    with open(restart) as f:
        ls=f.readlines()
    with open(restart,'w') as f:
        for l in ls:
            par=ls.strip('\n').split(':').strip(' ')
            if par[0] == 'totStep':
                stps=int(par[1])
                newsteps=stps+extsteps
                newl=par[0]+':'+str(newsteps)
                f.write(newl+'\n')
            else:
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
    os.system(exec_name+arg_start + ' > screen') 
else:#restart extending if requested, i.e. if totsteps > 0 
    if totsteps > 0:
        extend_sim(restart)
    os.system(exec_name + arg_restart + ' >> screen')
quit()
