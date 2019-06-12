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
# i nomi degli script di start e restart
import sys
import os
import psutil
from operator import itemgetter
#questo rescheduler è abbastanza portabile infatti funziona 
#sia in linux che in mac osx
#
def is_integer(s):
    try:
        int(s)
        return True
    except (ValueError, TypeError):
        return False
def print_error():
    print('You have to supply a config file where:')
    print('first lines=<tosteps=-1|>0> <extra steps=-1|>0> <max jobs>')
    print('list of jobs to check (with absolute paths)')
    print('reschedluer [-sv|-f<filter string>|-s/-show|-v/-verbos|-t/-type <type>|-extargs|-ea|-k/-kill] <conf file>')
    quit()
show_only = False
verbose = False   
args=sys.argv
sched_type='simpp'
filter_proc=''
extra_args=''
#print('args=',args)
if len(args)==1:
    print_error()
    quit()
lof=''
del(args[0])
killp=False
for a in args:
    if a == '-show' or a  == '-s':
        show_only=True
    elif a == '-verbose' or a == '-v':
        verbose=True    
    elif a == '-sv':    
        show_only=True
        verbose=True
    elif a == '-filter' or a == '-f':
        filter_proc=next(a)
    elif a == '-extargs' or a == '-ea':
        extra_args=next(a)
    elif a == '-type' or a == '-t':
        sched_type=next(a)
    elif a == '-kill' or a == '-k':
        killp=True
    else:
      lof = a
if lof == '':
    print_error()
    quit()
#print ('lof=',lof)    
totsteps=-1
extsteps=-1
#print ('uid=', os.getuid())
def can_access_cwd(pr):
    try:
        pr.cwd()
        return True
    except psutil.AccessDenied:
       return False 
def get_proc_info(fil):
    allpids=[]
    allcwds=[]
    allcls=[]
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
        cli=' '.join(pr.cmdline())
        if not can_access_cwd(pr):
            #print('qui cl=', pr.cmdline())
            continue
  #           include_dir=False
  #      else:
  #          include_dir=True    
        #print('pid=', pid, ' cl=', pr.cmdline())
        if len(fil)==0 or cli.find(fil)!=-1:
            #print('cli1=', cli, ' boh=', cli[1].find(fil))
            #print('pr.cwd=', pr.cwd())
            #print('proc pid=', pr.pid)
            allcls.append(cli)#command line con cui è stato eseguito			
            allpids.append(pr.pid)#pid del processo	
        #    if include_dir==True:
            allcwds.append(pr.cwd())#directory del processo
         #   else:
          #      allcwds.append('none')
            #print ('pid= ', pid, ' cl=', cli)
        #if pid == 89288 or pid==89367:
        #    print('pid=', pid, ' cwd=', pr.cwd(), ' cli=', cli)
    return (allcls,allpids,allcwds)
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
    return int(l[3])
#######################################
# VARIABILI CHE DIPENDONO DA TIPO DI PROGRAMMA DA RIAVVIARE
# SI PRESUPPONE COMUNQUE CHE ESISTANO DUE RESTART FILES
# CHE TERMINANO CON 0 o 1
#nomi dei restart files da valutare
postpend=' &'
prepend=''
#postpend=''
if sched_type == 'simpp':
    #restart can be also a list of one element!
    restart=['restart-0','restart-1']
    donefile='cnf-final' # se esiste questo file vuol dire che ha finito!
    #argomenti per l'eseguibile in caso
    # di start o restart
    arg_start=' 2'
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
# keep_going=False #if true restart finished simulations too
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
            #print('both restart files do not exist...I give up!')
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
def job_is_running(bn,en,clines,cwds):
    cc=0    
    for l in cwds:
#        if clines[cc].find(en)!=-1 and l == 'none':
#            print('Executable name ' + en + ' matches') 
#            print('but folder is not accessible,')
#            print('if it has been launched through mosrun command, please')
#            print('use a script to start the job')
#            return True
        if bn == l and clines[cc].find(en)!=-1:
            #print ('bn=', bn, ' cline=', clines[cc], ' en=', en)
            return True
        cc += 1
    return False
def to_extend(dir,w):        
    return totsteps > 0 and os.path.exists(dir+'/'+donefile) and w != -1
#
with open(lof) as f:
    lines=f.readlines() 
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
#print('cls=',allcls)
#print('fil=',filter_proc)
#print('allcwds=',allcwds)
#print('pids=',pids)
ll=lines[0].strip('\n').split(' ')
#print ('ll=',ll)
totsteps=int(ll[0])
extsteps=int(ll[1])
max_jobs=int(ll[2])
del(lines[0])
if killp == True:
    for l in allcwds:
        #print ('bn=',bnc)
        #print ('en=',enc)
        found=False
        for ll in lines:
            bn, en=os.path.split(ll.strip('\n'))
            #print ('bn=', bn, ' bnc=', l)
            if bn==l and allcls[nline].find(en) != -1: 
                found=True
        if found:
            print('nline=', nline, ' killing process ', pids[nline])
            os.system(' kill '+ str(pids[nline]))
        nline += 1    
    quit()        
for l in lines:
    bn, en=os.path.split(l.strip('\n'))
    if not job_is_running(bn,en,allcls,allcwds):
    #print ('job ', i, ' is missing')
        dir=bn
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
                exec=' ./'+en+build_arg_restart(nline,extra_args,which)
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
        exec=' ./'+en+build_arg_restart(nline,extra_args,which)
        s2e=prepend + exec + postpend 
        print('[extend simulation, dir=',bn,'] exec: ', s2e)
        os.system(s2e) 
        nj += 1
        lstextended.append(l[1])
        nrun +=1
#       
# stampare liste di running e riavviati
if len(lstrun) > 0 and show_only == False:
    #for l in lstrun:
    print ('list run=', lstrun)
#if len(lstdead) > 0:
#    print ('list dead=', lstdead)
#if len(lststart) > 0:
#    print ('list start=', lststart)
#if len(lstext) > 0:
#    print ('list extend=', lstext)
if show_only == True:
    if verbose == True:
        print ('R= running; F=finished')
    if verbose == True:
        for l in lines:
            if l in lstrun:
                print ('R ', end='')
            elif l in lstdone:
                print ('F ', end='')
            else:
                print ('  ', end='')
            print (l,end='')
    print ('tot running ='+str(len(lstrun))+'; tot finished=', len(lstdone))
    quit()
if not ok:
	print('Some jobs (#'+str(ndead)+') were dead and I had to restart them!')
if ndone < len(lines):
    #print ('ndone=',ndone, ' nrun=', nrun, ' ndead=', ndead)
    print('Jobs now running='+str(nrun),end='')
    print(', completed='+str(ndone)+'/'+str(len(lines))+' (max:'+str(max_jobs)+')')
else:
    print('All done here!')
