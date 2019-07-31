# la seguente funzione estrae dal Makefile le definizioni di macro che vengono assegnate alla variabile MACROS 
# e che vengono usate per la compilazione
import os.path as p
DIR_OF_THIS_SCRIPT = p.abspath( p.dirname( __file__ ) )    
def getmacros():
    with open(DIR_OF_THIS_SCRIPT+"/Makefile") as f:
        lines=f.readlines()
    fine=False
    macros_n=''
    for line in lines:
        llst=line.strip('\n').split(' ')
        if len(llst) > 1 and llst[0] == 'export' and (llst[1].find('MACROS')!=-1):
            #print ('llst=', llst)
            if len(llst) > 1:
                mac=llst[1].strip('\n').split('=') 
                #print ('mac=', mac, 'line[1]=', line[1])
                if len(mac) > 1:
                    llst2 = mac[1].split(')')
                    #print ('llst2=', llst2)
                    macros_n=llst2[0].replace('(','').replace('$','')
                    fine=True
                    break
        if fine==True:
            break
    macros=''
    for line in lines:
        llst = line.split('=')
        if llst[0] == macros_n and len(llst) > 1:
            macros = llst[1].strip('\n')
            break
    mactmp=macros.split('#')
    macros=mactmp[0]
    #with open('macros2.output','w') as f:
    #    f.write(macros)
   
    #aggiungo uno spazio prima e dopo altrimenti ci sono problemi
    alst=macros.split(' ')
    maclst = []
    for l in alst:
        l.replace(' ', '')
        if len(l) == 0:
            continue
        if l[1]=='U':
            maclst.append('-U')
        else:
            maclst.append('-D')
        lasti=len(l)
        maclst.append(l[2:lasti])
    #with open('macros.output','w') as f:
    #    f.write(str(maclst))
    return maclst

def Settings( **kwargs ):
    #language = kwargs[ 'language' ]
    #if language=='cfamily':
    #if language == 'python':
    #
    MACROS=getmacros()
    flags = [ '-x', 'c', '-Wall', '-Wextra', '-I', DIR_OF_THIS_SCRIPT+'/commSrc']
    #flags.append(' '+ MACROS +' ')
    flags=flags+MACROS
    return {'flags': flags}
