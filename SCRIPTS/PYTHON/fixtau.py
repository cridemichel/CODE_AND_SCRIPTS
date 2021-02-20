import os, sys
import numpy as np
args=sys.argv
if len(args) > 1:
    fn=args[1]
else:
    os.system('ls NPT*/triatpars.xasc NTV*/triatpars.xasc > _lista_')
    fn='_lista_'
print('filename=',fn)
K=np.sqrt(2.0/(2.0*np.pi)/(2.0*np.pi))
print('K=', K)
with open(fn) as f:
    fns=f.readlines()
lnew=[]
for pf in fns:
    os.system('cp '+pf.strip('\n')+ ' ' + pf.strip('\n') + '.bak' )
    with open(pf.strip('\n')) as f:
        lines=f.readlines()
        for l in lines:
            if l.find('#')!=-1:
                com=' #'+l.strip('\n').split('#')[1]
                lc=l.strip('\n').split('#')[0]
            else:
                com=''
                lc=l.strip('\n')
            ls=lc.split(':')
            if ls[0]=='taup' or ls[0]=='taub':
                lnew.append(ls[0]+': '+str(float(ls[1])*K)+com+'\n')
            else:
                lnew.append(l)
    with open(pf.strip('\n'),'w') as f:
        f.writelines(lnew)
