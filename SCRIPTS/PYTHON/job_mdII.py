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
prepend=' ' # prepend to exec name
exec_name=prepend+'nohup mosrun ./hard_cylinders ' # sim executable
#exec_name=prepend+'./hard_cylinders ' # sim executable
arg_start=' -fa ellipsoid_flex_mc.par ' # command line arg for starting
arg_restart=' -ca ' + restart #command line arg for restarting
##########################################################
# SIMULATION DEPENDENT CODE
# following functions depend on simulation file format
#
# get current simulation steps from restart file
def get_steps(restart):
    with open(restart) as f:
        ls=f.readlines()
    for l in ls:
        ll= l.strip('\n').split(':')
        if ll[0].strip().find('curStep')!=-1:
            curstps=int(ll[1])        
    return curstps
#
#extend simulation by modifying restart file
def extend_sim(restart):   
    #os.system('rm '+ donefile)
    with open(restart) as f:
        ls=f.readlines()
    cc=0
    for l in ls:
        ll= l.strip('\n').split(':')
        if ll[0].strip().find('totStep')!=-1:
            ls[cc] = ll[0]+': '+ str(int(ll[1])+extsteps) + '\n'
        cc+=1     
    #overwrite file total steps with new value
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
    os.system(exec_name+arg_start + ' > screen 2>&1') 
else:#restart extending if requested, i.e. if totsteps > 0 
    if totsteps > 0:
        extend_sim(restart)
    os.system(exec_name + arg_restart + ' >> screen 2>&1')
quit()
