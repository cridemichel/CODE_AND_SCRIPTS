#!/usr/bin/env python3
#questo script python assume che per ogni cartella fornita nella lista
#giri un solo processo
import sys
import os
import psutil
from operator import itemgetter
#questo rescheduler è abbastanza portabile infatti funziona 
#sia in linux che in mac osx
#
args=sys.argv
if len(args) > 1:
    lof=args[1]
else:
    print('You have to supply a file with all jobs to check (with absolute paths)')
    quit()
if len(args) > 2:
    filter_proc=args[2] #common pattern in exec filenames
else:
    filter_proc=''
if len(args) > 3:
    extra_args=args[3] #extra args for launching the executable
else:
    extra_args=''
if len(args) > 4:
    sched_type=args[4] # simpp or veff for now
else:
    sched_type='simpp' # simpp or veff for now
#
def get_proc_info(fil):
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
        cli=pr.cmdline()
        if len(cli) > 1:
            if fil == '' or cli[1].find(fil):
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
def get_steps(bn,w):
    with open(bn+'/'+restart[w]) as f:
        firstline=f.readline()
    l=firstline.split(' ')
    #print ('l=',l)
    #print ('lines[0][2]=', ls[0][2])
    #overwrite filed total steps with new value
    return int(l[2])
#######################################
# VARIABILI CHE DIPENDONO DA TIPO DI PROGRAMMA DA RIAVVIARE
# SI PRESUPPONE COMUNQUE CHE ESISTANO DUE RESTART FILES
# CHE TERMINANO CON 0 o 1
#nomi dei restart files da valutare
#prepend='/usr/bin/nohup /bin/mosrun '# li mette prima dell'exec
postpend=' >> screen &'#li mette dopo l'exec
prepend=''
#postpend=''
max_jobs=2
if sched_type == 'simpp':
    keep_going=True #if true restart finished simulations too
    #restart can be also a list of one element!
    restart=['restart-0','restart-1']
    donefile='cnf-final' # se esiste questo file vuol dire che ha finito!
    #argomenti per l'eseguibile in caso
    # di start o restart
    arg_start=' 2 '
    arg_restart=[' 1 restart-0 ', ' 1 restart-1 ']
    #maximum number of running jobs     
    #
    # NOTE SUL RESTART
    # se le simulazioni arrivano a totsteps i file di restart vengono cancellati
    # e non si riavviano più finendo senza riprendere
    # nel file di parametri delle sim vanno impostati gli step a cui termina 
    # la simulazione (steps) senza cancellare i restart e quelli a cui deve finire 
    # cancellando i restart per non proseguire.
    # Notare che quando raggiunge totsteps comunque scrive dei file di restart
    # chiamati restart-final-[0,1] che possono essere utilizzati per ripartire 
    # ulteriormente rinominandlo in restart-[0,1] 
    # NOTES:
    # customize these functions is special args, 
    # l is an integer between 0 and the total number of jobs
    # which is equal to len(lines) (=# of lines in the file supplied 
    # as argument)
    def build_arg_start(l,ea):
        return arg_start
    def build_arg_restart(l,ea,w):
        return arg_restart[w]
    # test if simulation is finished (customize if needed)
    def sim_done(dir):
        return os.path.exists(dir+'/'+donefile)        
    #step to extend simulation
    extra_steps=1000000
    #
    def extend_sim(bn,w):   
        os.system('rm '+donefile)
        with open(bn+'/'+restart[w]) as f:
            ls=f.readlines()
        l=ls[0].split(' ')
        #print ('l=',l)
        #print ('lines[0][2]=', ls[0][2])
        #overwrite filed total steps with new value
        newsteps=int(l[2])+extra_steps
        l[2]=str(newsteps)
        ls[0]=' '.join(l)
        #print ('dopo ls[0]=',ls[0])
        #print ('dir=', bn)
        #print ('res='+restart[w])
        with open(bn+'/'+restart[w],"w") as f:
            for l in ls:
                f.write(l)
elif sched_type == 'veff':
#la seguente va bene se il criterio è se ci sia un file
#con una certa linea finale 
#routine riavviare programmi con un solo file di restart 
# come quelli usati per il calcolo del potenziale efficace
# queste devono rimpiazzare le corrispondenti build_arg_start,
# build_arg_restart e sim_done
# inoltre anche la variabile restart_veff deve essere rinominata in restart
    keep_going=False #if true restart finished simulations too
    def sim_done(dir):
        with open(dir+'/veff_vs_tt.dat') as f:
            ls=f.readlines()
            lastline=ls[-1]
            lst=lastline.strip('\n').split(' ')
            if lst[0] == '99900000000':
                return True
            else:
                return False
            #            
            #build arg depending on l
    def build_arg_restart(l,ea,w):
        return '100000000000 ' + ea + ' 300 ' +str(l) + ' 100000000'
    def build_arg_start(l,ea):
        return '100000000000 ' + ea + ' 300 ' +str(l) + ' 100000000'
    restart=['calcveff.chk']
else:
    print('wrong schedule type '+sched_type)
    quit()
#######################################
# se il programma di restart è solo uno la seguente funzione 
# va cambiata opportunamente
def choose_restart(bn):
    if len(restart) == 2:
        f0n=bn+'/'+restart[0]
        f1n=bn+'/'+restart[1]
        ex0=os.path.exists(f0n)
        ex1=os.path.exists(f1n)
        #print("ex=", ex0, ex1)
        if ex0 == False and ex1 == False:
            print('both restart files do not exist...I give up!')
            return -1
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
        return which
    else:#only one restart file
        f0n=bn+'/'+restart[0]
        ex0=os.path.exists(f0n)
        if ex0 == False:
            return -1
        else:
            return 0            
#
with open(lof) as f:
    lines=f.readlines() 
c=0
#return command lines, pids and absolute path
cls,pids,allcwds=get_proc_info(filter_proc)
ok=True
ndone=0
ndead=0
nrun=0
#we compare absolute paths to determine if process is running
#(we assume that each jobs has been launched from a different dir)
nline=0
lstdone=[]
lststps=[]
#print('allcwds=',allcwds)
#print('pids=',pids)
for l in lines:
    bn, en=os.path.split(l.strip('\n'))
    if bn not in allcwds:
    #print ('job ', i, ' is missing')
        dir=bn
        if sim_done(dir): 
            if keep_going == True:
                which=choose_restart(bn)
                lstdone.append([get_steps(dir,which),l])
            ndone+=1 
        else:
            ndead+=1
            if ndead + nrun > max_jobs:
                print ('maximum number of jobs (' + str(max_jobs) + ') reached')
                continue
            print('job '+ en + ' is not running and it has not finished yet!', end='')
            print(' I am restarting it...')
            which=choose_restart(bn)
            if which == -1:
                #no restart file, first start
                exec=' ./'+en+build_arg_start(nline,extra_args)
                os.chdir(dir)
                #print('dir=',os.getcwd())
                s2e=prepend + exec + postpend 
                print('executing ', s2e)
                os.system(s2e)
            else:
                #print('en=',en)
                exec=' ./'+en+build_arg_restart(nline,extra_args,which)
                print ('exec is: '+exec)
                os.chdir(dir)
                #print('dir=',os.getcwd())
                s2e=prepend + exec + postpend 
                print('executing ', s2e)
                os.system(s2e)
                ok=False
    else:#bn is running if here
        nrun+=1
    nline+=1
if keep_going == True:
    #sort list of done simulations from lower steps to higher
    #so that simulations left behind start first
    lstdone.sort(key=itemgetter(0))
    nj=0
    njmax = max_jobs-(nrun+ndead)
    for l in lstdone:
        if nj >= njmax:
            break
        bn, en=os.path.split(l[1].strip('\n'))
        os.chdir(bn)
        which=choose_restart(bn)
        extend_sim(bn,which)
        exec=' ./'+en+build_arg_restart(nline,extra_args,which)
        s2e=prepend + exec + postpend 
        print('[keepgoing] executing ', s2e)
        os.system(s2e) 
        nj += 1
#       
if not ok:
	print('Some jobs (#'+str(ndead)+') were dead and I had to restart them!\n')
else:
    if keep_going == True:
        print ('nj=',nj, ' nrun=', nrun, ' pndead=', ndead)
        print('[keepgoing] Jobs runing='+str(nj+nrun+ndead))
        print('max_jobs='+str(max_jobs)+ ' total jobs='+str(len(lines)))
    else:
        if ndead == 0 and ndone == len(lines):
            print('All done here!')
        else:
            print('There are '+str(nrun)+' jobs runnings', end='')
            print(' and '+ str(ndone) + ' regularly finished')
            print('Total number of job is', str(len(lines)))
