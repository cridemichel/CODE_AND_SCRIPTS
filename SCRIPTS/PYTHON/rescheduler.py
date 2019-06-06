#!/usr/bin/env python3
import sys
import os
import psutil
#questo rescheduler è abbastanza portabile infatti funziona 
#sia in linux che in mac osx
def get_proc_cmdlines():
    allpids=[]
    allcwds=[]
    cls=[]
    pids = [pid for pid in psutil.pids()]
    for pid in pids:
        pr=psutil.Process(pid)
        uid=pr.uids()[1]#effective uid(0 is uid, 1 is effective uid)
        #lo uid potrebbe essere quello dell'utente mentre l'effettivo
        #potrebbe essere 0 se si è usato setuid cosicché 
        #si ottiene un permission denied poiché il processo è rooted
        #quindi va usato l'effective uid
        if uid != os.getuid():
        #considera solo i processi che appartengono all'utente
            continue
        #print ('p=',p, 'cmd=',pr.cmdline())
        #filtra le command line con la stringa 
        cls.append(pr.cmdline())#command line con cui è stato eseguito			
        allpids.append(pr.pid)#pid del processo	
        allcwds.append(pr.cwd())#directory del processo
    return (cls,allpids,allcwds)
def get_num_words(fn):
    with open(fn) as f:
        lines=f.readlines()
        nw=0
        for line in lines:
                l=line.strip('\n').split(' ')
                nw+=len(l)
        return nw
#######################################
# VARIABILI CHE DIPENDONO DA TIPO DI PROGRAMMA DA RIAVVIARE
# SI PRESUPPONE COMUNQUE CHE ESISTANO DUE RESTART FILES
# CHE TERMINANO CON 0 o 1
#nomi dei restart files da valutare
prepend='/usr/bin/nohup /bin/mosrun '# li mette prima dell'exec
postpend=' >> screen &'#li mette dopo l'exec
restart0='restart-0'
restart1='restart-1'
donefile='cnf-final' # se esiste questo file vuol dire che ha finito!
#######################################
arg=' 1 restart-' #argomenti per l'eseguibile
args=sys.argv
if len(args) > 1:
    lof=args[1]
else:
    print('You have to supply a file with all jobs to check')
    quit()
#if len(args) > 2:
#    cfn=args[2] #common patter in exec filenames
#else:
#    cfn=''
with open(lof) as f:
    lines=f.readlines() 
c=0
#return command lines, pids and absolute path
cls,pids,allcwds=get_proc_cmdlines()
ok=True
ndone=0
ndead=0
nrun=0
#we compare absolute paths to determine if process is running
#(we assume that each jobs has been launched from a different dir)
for l in lines:
    bn, en=os.path.split(l.strip('\n'))
    if bn not in allcwds:
    #print ('job ', i, ' is missing')
        dir=bn
        if os.path.exists(dir+'/'+donefile):
            ndone+=1 
        else:
            ndead+=1
            print('job '+ en + ' is not running and it has not finished yet!', end='')
            print(' I am restarting it...')
            f0n=bn+'/'+restart0
            f1n=bn+'/'+restart1
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
            exec=' ./'+en+arg+str(which)
            print ('exec is: '+exec)
            os.chdir(dir)
            print('dir=',os.getcwd())
            s2e=prepend + exec + postpend 
            print('executing ', s2e)
            os.system(s2e)
            ok=False
    else:#bn is running if here
        nrun+=1
if not ok:
	print('Some jobs (#'+str(ndead)+') were dead and I had to restart them!\n')
else:
	if ndead == 0 and ndone == len(lines):
		print('All done here!')
	else:
		print('There are '+str(nrun)+' jobs runnings', end='')
		print(' and '+ str(ndone) + ' regularly finished')
		print('Total number of job is', str(ndone+nrun))
