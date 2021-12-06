#!/usr/bin/env python3
#questo script python assume che per ogni cartella fornita nella lista
#giri un solo processo
# N.B. mosrun viene eseguito con setuid root quindi le info del processo
# (ad es. il path) non sono accessibili per cui se si usa mosrun 
# è bene lanciare la simulazione tramite uno script che poi lancia l'eseguibile con mosrun
# Per passare gli argomenti basta usare $@ per passare tutti gli arg dello script all'eseguibile
# lanciato da mosrun.
# Ora questo script python presuppone l'esistenza di un'unico script/eseguibuile
# che viene usato per lanciare e rilanciare, si potrebbe anche pensare di creare
# un file di configurazione solo con le directory e di mettere nella prima riga
# i nomi degli script di start e restart.
import sys
import os
import psutil
import signal
from operator import itemgetter
found_itutil=True
try:
    from iteration_utilities import unique_everseen
    from iteration_utilities import duplicates
except ImportError:
    print('Module iteration_utilities not found')
    print('You can install it with:')
    print('> pip3 install iteration_utilities')
    print('I am not gonna use it now...')
    found_itutil=False
#questo rescheduler è abbastanza portabile infatti funziona 
#sia in linux che in mac osx
#
#initializations do not touch!
restart=[]
donefile=''
jobfinished=''
def is_integer(s):
    try:
        int(s)
        return True
    except (ValueError, TypeError):
        return False
def print_error():
    print('rescheduler [-sv|-filter <filter string>|-s/-show|-v/-verbos|-t/-type <type>|-extargs|-ea|-k/-kill <list_to_kill>|')
    print('-delete|-d <list_to_delete>|-finish|-f <list_to_finish>|-kf <list_to_kill_and_finish>')
    print('|-kd <list_to_kill_and_delete>|-kk <conf_file>')
    print('where <list_to_kill>=0,1,2,3-4 (if equal to \'all\' means all jobs)')
    print('-kf kill specified processes and mark them as done')
    print('-kd kill specified processes and delete job from configuration file')
    print('-kk if a process does not terminate with SIGTERM use SIGKILL signal')
    print('and <conf_file> is a configuration file with the following structure:\n')
    print('<tosteps=-1|>0> <extra steps=-1|>0> <max jobs> <donefile> <jobfinished> <restart0> <restart1>')
    print('/home/demichel/jobs1.sh\n/home/demichel/jobs2.sh')
    print('\n<totsteps> is the total number of steps (-1 means to not extend)')
    print('<extra steps> is the number of steps to extend simulations')
    print('<donefile> is the file written when a simulation ends (if totsteps < 0 jobs is finished)')
    print('<jobfinished> is the file written when jobs is terminated and, if it exists, job will not be restarted anymore')
    print('<restart0> and <restart1> are the names of the restart files')
    print('<donefile>, <jobsfinished>, <restart0> and <restart1> are optional')
    print('first line is followed by a list of jobs (preferably shell scripts) to check with absolute paths')
    quit()
list_to_kill=[int]
list_to_finish=[int]
list_to_delete=[int]
def parse_ranges(arg):
    lst=[]
    sarg=[]
    if arg != 'all':
        if arg.find(',')!=-1:
            sarg=arg.split(',')
        else:
            sarg.append(arg)
    else: 
        return []
    for l in sarg:
        ll=l.split('-')
        #print('ll=',ll, ' arg=' , arg)
        if len(ll)==1:
            if not ll[0].isnumeric():
                print_error()
                quit()
            lst.append(int(ll[0]))
        elif len(ll)==2:
            if (not ll[0].isnumeric()) or (not ll[1].isnumeric()):
                print_error()
                quit()
            for i in range(int(ll[0]),int(ll[1])+1):
                lst.append(i)
            #print('lst=', lst)
        else:
            print('Parse error in list of jobs to kill')
            quit()
    return lst
def check_range(lista,lines):
    for i in lista:
        if i < 0 or i >= len(lines):
            print('Job index out of range [0,'+ str(len(lines)-2) +']')
            quit()
if found_itutil:
    def check_duplicates(lstchk):
        retlst=list(unique_everseen(duplicates(lstchk)))
        if len(retlst) > 0:
            print('The following entries are duplicated:')
            for lll in retlst:
                print(lll.strip('\n'))
            print('Please fix it!')
            quit()

show_only = False
verbose = False
args=sys.argv
sched_type=''
filter_proc=''
extra_args=''
#print('args=',args)
if len(args)==1:
    print_error()
    quit()
lof=''
del(args[0])
killp=False
karg=''
itargs=iter(args)
deljobs=False
finjobs=False
restartjobs=False
usesigkill=False
for a in itargs:
    if a == '-show' or a  == '-s':
        show_only=True
    elif a == '-verbose' or a == '-v':
        verbose=True
    elif a == '-kk':
        usesigkill=True
    elif a == '-sv':
        show_only=True
        verbose=True
    elif a == '-filter':
        filter_proc=next(itargs)
    elif a == '-extargs' or a == '-ea':
        extra_args=next(itargs)
    elif a == '-type' or a == '-t':
        sched_type=next(itargs)
    elif a == '-kill' or a == '-k':
        karg=next(itargs)
        killp=True
    elif a == '-delete' or a == '-d':
        darg=next(itargs)
        deljobs=True
    elif a == '-finish' or a == '-f':
        farg=next(itargs)
        finjobs=True
    elif a == '-restart' or a == '-r':
        rarg=next(itargs)
        restartjobs=True
    elif a== '-kf':
        karg=next(itargs)
        killp=True
        finjobs=True
    elif a=='-kd':
        karg=next(itargs)
        killp=True
        deljobs=True
    else:
        lof = a
if not os.path.exists(lof):
    print('file '+ lof + ' does not exist')
    quit()
with open(lof) as f:
    lines=f.readlines()
if found_itutil:
    check_duplicates(lines)
if killp==True:
    list_to_kill=parse_ranges(karg)
    check_range(list_to_kill,lines)
    #print('karg=',karg)
    #print("list_to_kill=", list_to_kill)
    #quit()        
if deljobs==True:
    if killp == True:
        list_to_delete=list_to_kill
    else:
        list_to_delete=parse_ranges(darg)
        check_range(list_to_delete,lines)
if restartjobs==True:
    list_to_restart=parse_ranges(rarg)    
if finjobs==True:
    if killp == True:
        list_to_finish=list_to_kill
    else:
        list_to_finish=parse_ranges(farg)
        check_range(list_to_finish,lines)
if lof == '':
    print_error()
    quit()
#print ('lof=',lof)    
totsteps=-1
extsteps=-1
#print ('ll=',ll)
ll=lines[0].strip('\n').split()
try:
    totsteps=int(ll[0])
    extsteps=int(ll[1])
    max_jobs=int(ll[2])
except ValueError:
    print('Value in config file is not an integer, exiting...')
    quit()
if (len(ll) > 3):
    donefile=ll[3]
if (len(ll) > 4):
    jobfinished=ll[4]
if (len(ll) > 5):
    restart.append(ll[5])
if (len(ll) > 6):
    restart.append(ll[6])  
firstline=lines[0]
del(lines[0])
#print ('uid=', os.getuid())
#return a tuple where first value is True if directory has been accessed
#and False otherwise while second value is an integer where 0 means
#that process died meanwhile can not be appended to the list
def can_access_cwd(pr):
    try:
        pr.cwd()
        return True
    except (psutil.AccessDenied,psutil.NoSuchProcess):
        return False
def get_proc_info(fil):
    allpids=[]
    allcwds=[]
    allcls=[]
    pids = [pid for pid in psutil.pids()]
    for pid in pids:
        try:
            pr=psutil.Process(pid)
        except psutil.NoSuchProcess:
            continue
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
        if not can_access_cwd(pr):
            #print('qui cl=', pr.cmdline())
            #          continue
            include_dir=False
        else:
            include_dir=True    
        try:
            cli=' '.join(pr.cmdline())
            prid = pr.pid
            if include_dir==True:
                pdir = pr.cwd()
            else:
                pdir='none'
        except (psutil.NoSuchProcess):
            continue
        #print('pid=', pid, ' cl=', pr.cmdline())
        if len(fil)==0 or cli.find(fil)!=-1:
            #print('cli1=', cli, ' boh=', cli[1].find(fil))
            #print('pr.cwd=', pr.cwd())
            #print('proc pid=', pr.pid)
            allcls.append(cli)#command line con cui è stato eseguito			
            allpids.append(prid)#pid del processo	
            allcwds.append(pdir)#directory del processo
            #print ('pid= ', pid, ' cl=', cli)
        #if pid == 89288 or pid==89367:
        #    print('pid=', pid, ' cwd=', pr.cwd(), ' cli=', cli)
    return (allcls,allpids,allcwds)
def get_num_words(fn):
    with open(fn) as f:
        lines=f.readlines()
        nw=0
        for line in lines:
            l=line.strip('\n').split()
            nw+=len(l)
        return nw
#######################################
# VARIABILI CHE DIPENDONO DA TIPO DI PROGRAMMA DA RIAVVIARE
# SI PRESUPPONE COMUNQUE CHE ESISTANO DUE RESTART FILES
# CHE TERMINANO CON 0 o 1
#nomi dei restart files da valutare
postpend=' &'
prepend=''
#postpend=''
if sched_type == 'simpp':
    def get_steps(bn,w):
        with open(bn+'/'+restart[w]) as f:
            firstline=f.readline()
            l=firstline.split()
        #print ('l=',l)
        #print ('lines[0][2]=', ls[0][2])
        #overwrite filed total steps with new value
        return int(l[3])
    #restart can be also a list of one element!
    restart=['restart-0','restart-1']
    donefile='cnf-final' # se esiste questo file vuol dire che ha finito!
    #argomenti per l'eseguibile in caso
    # di start o restart
    arg_start=' 2'
    arg_restart=[' 1 restart-0 ', ' 1 restart-1 ']
    jobfinished='_DONE_'
    # maximum number of running jobs     
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
    def build_arg_restart(l,ea,w,ext):
        return arg_restart[w]
    # test if simulation is finished (customize if needed)
    def sim_done(dir):
        if os.path.exists(dir+'/'+jobfinished):
            return True
        if totsteps > 0:
            w=choose_restart(dir)
            s=get_steps(dir,w)
            if s >= totsteps:
                return True
            else:
                return False
        else:     
            return os.path.exists(dir+'/'+donefile)        
    #step to extend simulation
    #extra_steps=1000000
    #
    def extend_sim(bn,w):   
        os.system('rm '+donefile)
        with open(bn+'/'+restart[w]) as f:
            ls=f.readlines()
        l=ls[0].split()
        #print ('l=',l)
        #print ('lines[0][2]=', ls[0][2])
        #overwrite filed total steps with new value
        newsteps=int(l[2])+extsteps
        l[2]=str(newsteps)
        ls[0]=' '.join(l)
        #print ('dopo ls[0]=',ls[0])
        #print ('dir=', bn)
        #print ('res='+restart[w])
        with open(bn+'/'+restart[w],"w") as f:
            for l in ls:
                f.write(l)
elif sched_type == 'veff':
    # la seguente va bene se il criterio è se ci sia un file
    # con una certa linea finale 
    # routine riavviare programmi con un solo file di restart un solo file di restart 
    # come quelli usati per il calcolo del potenziale efficace
    # queste devono rimpiazzare le corrispondenti build_arg_start,
    # build_arg_restart e sim_done
    jobfinished='_DONE_'
    def get_steps(bn,w):
        return 0 
    restart=['calcveff.chk']
    donefile='not_used' # se esiste questo file vuol dire che ha finito!
    #argomenti per l'eseguibile in caso
    # di start o restart
    arg_start='not_used'
    arg_restart=['not_used', 'not_used']
    def sim_done(dir):
        if os.path.exists(dir+'/'+jobfinished):
            return True
        with open(dir+'/veff_vs_tt.dat') as f:
            ls=f.readlines()
            lastline=ls[-1]
            lst=lastline.strip('\n').split()
            if lst[0] == '99900000000':
                return True
            else:
                return False
            #            
            #build arg depending on l
    def build_arg_restart(l,ea,w,ext):
        return '100000000000 ' + ea + ' 300 ' +str(l) + ' 100000000'
    def build_arg_start(l,ea):
        return '100000000000 ' + ea + ' 300 ' +str(l) + ' 100000000'
    def extend_sim(bn,w):
        pass
else:
    if len(restart)==0:
        restart=['restart-0','restart-1']
    # donefile is the one created when a run terminates 
    if donefile=='':
        donefile='cnf-final'
    #get modificatiom time (which is a float) instead of steps in this general case
    def get_steps(bn,w):
        #print('mod time=',os.path.getmtime(os.path.join(dir,donefile)))
        return os.path.getmtime(os.path.join(dir,donefile))
    #
    if jobfinished=='':
        jobfinished='_DONE_'
    # test if simulation is finished (customize if needed)
    # [totsteps > 0] if the scripts creates a file called DONE simulation 
    #   will not be restarted 
    # [totsteps <= 0] if a file called donefile is created jobs is finished
    #print ('jobfib=',jobfinished)
    def sim_done(dir):
        if os.path.exists(dir+'/'+jobfinished):
            return True
        if totsteps <= 0:
            return os.path.exists(dir+'/'+donefile)
        return False
    #ea are extra arg which can be specified in the command line
    # if ext=True provide totsteps and exsteps otherwise pass -1 -1 to script 
    def build_arg_restart(l,ea,w,ext):
        if ext==True:
            return ' 1 '+restart[w]+' '+str(totsteps)+' '+str(extsteps)+' '+jobfinished+' '+ea
        else:
            return ' 1 '+restart[w]+' -1 -1 '+jobfinished+' '+ea
    def build_arg_start(l,ea):
        return ' 2 '+ restart[0] + ' -1 -1 ' + jobfinished + ' ' + ea #2 means first start 
    def extend_sim(bn,w):
        pass #not used
    #print('restart=', restart, ' donefile=', donefile)
    #print('wrong schedule type '+sched_type)
    #quit()
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
            #print('both restart files do not exist...I give up!')
            return -1
        if ex0 == False:
            which=1
        elif ex1 == False:
            which=0
        else:
            nw0=get_num_words(f0n)
            nw1=get_num_words(f1n)
            #print('nw0=', nw0, ' nw1=', nw1)
            if nw0 > nw1:
                which=0
            elif nw0 < nw1:
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
def job_is_running(bn,en,clines,cwds):
    cc=0    
    for l in cwds:
        if clines[cc].find(en)!=-1 and l == 'none':
            print('Executable name ' + en + ' matches') 
            print('but folder is not accessible,')
            print('if it has been launched through mosrun command, please')
            print('use a script to start the job and change ')
            print('the configuration file accordingly')
            return True
        if bn == l and clines[cc].find(en)!=-1:
            #print ('bn=', bn, ' cline=', clines[cc], ' en=', en)
            return True
        cc += 1
    return False
def to_extend(dir,w):        
    return totsteps > 0 and os.path.exists(dir+'/'+donefile) and w != -1
#
def on_terminate(proc):
    print("process {} terminated with exit code {}".format(proc, proc.returncode))

def kill_parent_and_childs(pid):
    parent = psutil.Process(pid)    
    children = parent.children(recursive=True)
    children.append(parent)#kill also the parent
    for p in children:
        try:
            p.send_signal(signal.SIGTERM)
        except (psutil.NoSuchProcess,psutil.AccessDenied):
            pass
    gone, alive = psutil.wait_procs(children,timeout=3,callback=on_terminate)
    if usesigkill==True and alive:
        for p in alive:
            print("process {} survived SIGTERM; trying SIGKILL".format(p))
            try:
                p.kill()
            except (psutil.NoSuchProcess,psutil.AccessDenied):
                pass
        gone, alive = psutil.wait_procs(alive, timeout=3, callback=on_terminate)
        if alive:
            # give up
            for p in alive:
                print("process {} survived SIGKILL; giving up".format(p))
c=0
#return command lines, pids and absolute path
allcls,pids,allcwds=get_proc_info(filter_proc)
#print ('allcwd=', allcwds)
#print ('allpids=', pids)
ok=True
ndone=0# numero terminati
ndead=0#numero morti
nrun=0#numero running
#we compare absolute paths to determine if process is running
#(we assume that each jobs has been launched from a different dir)
nline=int(0)
lstdone=[]
lstrun=[]
lstdead=[]
lststart=[]
lsttoext=[]
lstextended=[]
lst_not_exist=[]
#print('cls=',allcls)
#print('fil=',filter_proc)
#print('allcwds=',allcwds)
#print('pids=',pids)
if killp == True:
    #pids = [pid for pid in psutil.pids()] 
    nkilled=0
    lst_kill_found=[]
    for l in allcwds:
        #print ('bn=',bnc)
        #print ('en=',enc)
        found=False
        cc=0
        #print('listtokill=', list_to_kill)
        for al in lines:
            if len(list_to_kill) > 0 and (not (cc in list_to_kill)):
                cc+=1
                continue
            bn, en=os.path.split(al.strip('\n'))
            #print ('bn=', bn, ' bnc=', l)
            if bn==l and allcls[nline].find(en) != -1: 
                found=True
                lst_kill_found.append(cc)
            if allcls[nline].find(en)!=-1 and l == 'none':
                print('Executable name ' + en + ' matches') 
                print('but folder is not accessible,')
                print('if it has been launched through mosrun command, please')
                print('use a script to start the job and change ')
                print('the configuration file accordingly')
                print('I skip it...')
                #found=False
            cc+=1
        if found:
            print('Killing process', pids[nline],end='')
            print(' and, recursively, all its subprocesses:')
            #os.system(' kill '+ str(pids[nline]))
            #kills all processes which have current pid as parent pid
            # suspend the job first
            try:
                #psutil.Process(pids[nline]).suspend()
                # kill its childs (so that if job is a script it kills all 
                # subprocesses)
                kill_parent_and_childs(pids[nline])
                # finally kill the job
                #psutil.Process(pids[nline]).terminate()
            except (psutil.NoSuchProcess,psutil.AccessDenied):
                continue
            nkilled+=1
        nline += 1    
    if nkilled == 0:
        print('No process running...')
    else:
        print (nkilled, 'processes killed together with all their childs')        
    if finjobs==False and deljobs==False: 
        quit()        
if restartjobs == True:
    cc=0 
    #print('list found=', lst_kill_found)
    #print ('qui')
    for l in lines:
        bn, en=os.path.split(l.strip('\n'))
        if not os.path.exists(bn):
            print('folder '+bn + ' does not exist, I skip it...')
            cc+=1
            continue
        if len(list_to_kill) == 0 or (cc in list_to_finish):
            os.chdir(bn)
            #if jobfinished files exists job is done
            if os.path.exists(jobfinished):
                os.system('rm '+ jobfinished)
        cc+=1 
    quit()
if finjobs == True:
    if killp == True:
        list_to_finish = list_to_kill
    cc=0 
    #print('list found=', lst_kill_found)
    #print ('qui')
    for l in lines:
        bn, en=os.path.split(l.strip('\n'))
        if not os.path.exists(bn):
            print('folder '+bn + ' does not exist, I skip it...')
            cc+=1
            continue
        if len(list_to_kill) == 0 or (cc in list_to_finish):
            if not os.path.exists(bn):
                print('folder '+bn + ' does not exist, I skip it...')
                cc+=1
                continue
            os.chdir(bn)
            #if jobfinished files exists job is done
            if not os.path.exists(jobfinished):
                os.system('touch '+jobfinished)
        cc+=1 
    quit()
if deljobs == True:
    if killp==True:
        list_to_delete = list_to_kill
    cc=0
    #print('ld=', list_to_delete)
    #backup
    os.system('cp ' + lof + ' ' + lof + '.bak')
    with open(lof,"w") as f: 
        f.write(firstline)
        for l in lines:
            if (not len(list_to_delete) == 0) and (not cc in list_to_delete):    
               f.write(l)
            cc+=1
    quit()
for l in lines:
    bn, en=os.path.split(l.strip('\n'))
    if not job_is_running(bn,en,allcls,allcwds):
        #print ('job ', i, ' is missing')
        dir=bn
        if not os.path.exists(dir):
            print('folder '+bn + ' does not exist, I skip it...')
            lst_not_exist.append(l)
            continue
        if sim_done(dir): 
            #job is finished correctly if here 
            ndone += 1 
            lstdone.append(l)
        else:
            if show_only == True:
                continue
            #not a cnf-final found (i.e. simulation is dead)
            #look for restart files
            if nrun >= max_jobs:
                #print ('maximum number of jobs (' + str(max_jobs) + ') reached')
                break
            which=choose_restart(bn)
            # se totsteps > 0 e ci sono il donefiled e il restart, prova ad estenderla
            # mettendola nella lista che poi verrà utilizzata per i riavvii
            if to_extend(dir,which):
                lsttoext.append([get_steps(dir,which),l])
                continue
            if which == -1:
                # se non ci sono i file di restart ma il donefile 
                # assumo che siano stati cancellati o rinominati per
                # per far terminare la simulazione
                if totsteps > 0 and os.path.exists(dir+'/'+donefile):
                    continue
                lststart.append(l)
                nrun += 1
                #no restart file, first start
                exec=' ./'+en+build_arg_start(nline,extra_args)
                os.chdir(dir)
                print('[first start, dir=',os.getcwd(),'] ',end='')
                s2e=prepend + exec + postpend 
                print('exec: ', s2e)
                os.system(s2e)
            else:
                ndead += 1
                lstdead.append(l)
                #at least one restart file found
                exec=' ./'+en+build_arg_restart(nline,extra_args,which,False)
                #print ('exec is: '+exec)
                os.chdir(dir)
                print('[restart dead, dir=',os.getcwd(),'] ',end='')
                s2e=prepend + exec + postpend 
                print('exec: ', s2e)
                nrun += 1
                os.system(s2e)
                ok=False
    else:#bn is running if here
        nrun+=1
        #lstrun.append([pids[nline],l])
        lstrun.append(l)
    nline+=1
# all jobs with a cnf-final are extended now 
# if number of jobs does not exceed max_jobs
if totsteps > 0 and show_only==False:
    #sort list of done simulations from lower steps to higher
    #so that simulations left behind start first
    lsttoext.sort(key=itemgetter(0))
    nj=0
    for l in lsttoext:
        if nrun >= max_jobs:
            break
        bn, en=os.path.split(l[1].strip('\n'))
        os.chdir(bn)
        which=choose_restart(bn)
        extend_sim(bn,which)
        exec=' ./'+en+build_arg_restart(nline,extra_args,which,True)
        s2e=prepend + exec + postpend 
        print('[extend simulation, dir=',bn,'] exec: ', s2e)
        os.system(s2e) 
        nj += 1
        lstextended.append(l[1])
        nrun +=1
#       
# stampare liste di running e riavviati
#if len(lstrun) > 0 and show_only == False:
    #for l in lstrun:
#    print ('list run=', lstrun)
#if len(lstdead) > 0:
#    print ('list dead=', lstdead)
#if len(lststart) > 0:
#    print ('list start=', lststart)
#if len(lstext) > 0:
#    print ('list extend=', lstext)
if show_only == True:
    if verbose == True:
        print ('R= running; F=finished; X=directory not found')
    if verbose == True:
        cc=0
        for l in lines:
            print ('[',str(cc),']',end='')
            if l in lstrun:
                print (' R ', end='')
            elif l in lstdone:
                print (' F ', end='')
            elif l in lst_not_exist:
                print (' X ', end='')
            else:
                print ('   ', end='')
            print (l,end='')
            cc+=1
    print ('tot running= '+str(len(lstrun))+'; tot finished=', len(lstdone),'/',len(lines), ' (max:'+str(max_jobs)+')')
    quit()
if not ok:
	print('Some jobs (#'+str(ndead)+') were dead and I had to restart them!')
if ndone < len(lines):
    #print ('ndone=',ndone, ' nrun=', nrun, ' ndead=', ndead)
    print('Jobs now running='+str(nrun),end='')
    print(', completed='+str(ndone)+'/'+str(len(lines))+' (max:'+str(max_jobs)+')')
else:
    print('All done here!')
