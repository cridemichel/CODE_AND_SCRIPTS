#!/usr/bin/env python3
import sys
import os
def get_proc_cmdlines(p=''):
    cls=[]
    allpids=[]
    pids = [pid for pid in os.listdir('/proc') if pid.isdigit()]
    for pid in pids:
        try:
            with open(os.path.join('/proc', pid, 'cmdline'), 'rb') as f:
                l=f.readline()
                l=l.replace(b'\x00',b'\x20')
                l2=str(l.decode('utf-8'))
                if p in l2:
                    cls.append(l2)			
                    allpids.append(pid)	
        except IOError: # proc has already terminated
            continue
    return (cls,allpids)
def get_num_words(fn):
    with open(fn) as f:
        lines=f.readlines()
        nw=0
        for line in lines:
                l=line.strip('\n').split(' ')
                nw+=len(l)
        return nw
#get X0
args=sys.argv
if len(args) > 1:
    lof=args[1]
else:
    print('You have to supply a file with all jobs to check')
    quit()
if len(args) > 2:
    cfn=args[2] #common patter in exec filenames
else:
    cfn=''
with open(lof) as f:
    lines=f.readlines() 
c=0
#return all pid obtained from commandlines containing string p
jobs,pids=get_proc_cmdlines(cfn)
#print("pids=",pids)
alljobs=[]
for e in jobs:
    l=e.split(' ')
		#print("process",pids[c],"=",l)
    alljobs.append(l)
#flatten the nested list making it a list from a nested list
alljobs=[item for sublist in alljobs for item in sublist]
#print ("all jobs=", alljobs)
ok=True
nd=0
for l in lines:
    bn, en=os.path.split(l.strip('\n'))
    rn='./'+en
    print ('rn=',rn)
    if rn not in alljobs:
    #print ('job ', i, ' is missing')
        dir=bn
        if os.path.exists(dir+'/cnf-final'):
            ok=True
        else:
            nd+=1
            print('job '+ en + ' is not running and it has not finished yet!', end='')
            print(' I am restarting it...')
            f0n=bn+'/restart-0'
            f1n=bn+'/restart-1'
            ex0=os.path.exists(f0n)
            ex1=os.path.exists(f1n)
            #print("ex=", ex0, ex1)
            if ex0 == False and ex1 == False:
                print('both restart files do not exist...I give up!')
                continue
            if ex0 == False: 
                which=1
            elif ex1 == False:
                which=0
            else: 
                nw0=get_num_words(f0n)
                nw1=get_num_words(f1n)
                if nw0 > nw1:
                    which=0
                elif nw1 < nw0:
                    which=1
                else:
                    f0 = os.path.getmtime(f0n)
                    f1 = os.path.getmtime(f1n)
                    if (f0 < f1):
                        which=1
                    else:
                        which=0
            #print('en=',en)
            exec=' ./'+en+' 1 restart-'+str(which)
            print ('exec is: '+exec)
            os.chdir(dir)
            print('dir=',os.getcwd())
            s2e=exec + ' >> screen &'
            print('executing ', s2e)
            os.system(s2e)
            ok=False
if not ok:
	print('Some jobs (#'+str(nd)+') were dead and I had to restart them!\n')
else:
	if c == 0:
		print('All done here!')
	else:
		print('There are '+str(c)+' jobs runnings', end='')
		print(' and '+ str(nd) + ' regularly finished')
		print('Total number of job is', str(nd+c))
