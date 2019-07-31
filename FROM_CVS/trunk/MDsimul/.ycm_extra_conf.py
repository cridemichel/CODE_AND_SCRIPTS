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
    with open('macros2.output','w') as f:
        f.write(macros)
   
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
    with open('macros.output','w') as f:
        f.write(str(maclst))
    return maclst

def Settings( **kwargs ):
    #language = kwargs[ 'language' ]
    #if language=='cfamily':
    #if language == 'python':
    #
    #MACROS=' -UMC_PERWER -DDIST_OF_CLOSE_APPR  -DPOLYELLIPS -UMC_HYDROPHOBIC_INT -DMC_NVE -UMC_FREEZE_BONDS -DMC_SOFTHE -UMC_SWELL -UMC_CLUSTER_MOVE -UMC_CLUSTER_NPT -UMC_ALT_ROT -DMC_FLIP_MOVE -UMC_OPT_BUILDATOMPOS -DMCIN_OPT -UMC_HC -DMC_SUS -UMC_CALC_COVADD -DMD_NO_SYSTEM -DMC_GRANDCAN  -DMC_SIMUL -UMD_STANDALONE -DMD_LXYZ -UMD_CHAIN_SIM -UMD_CALC_VBONDING -UDMD_MULTIPLE_LL -DMD_DYNAMIC_OPROG -UMD_CALENDAR_HYBRID -DMD_SCALEPHI_STAGES -DMD_OPT_SCALEPHI -UMD_ALLOW_ONE_IGG_BOND -UMD_PROTEIN_DESIGN -UMD_PARANOID_CHECKS -DMD_GRAZING_TRYHARDER -DMD_GHOST_IGG -DMD_RABBIT -DMD_SAVE_SPOTS -DMD_SUPERELLIPSOID -UMD_ABSORPTION -UMD_FOUR_BEADS -DMD_BASIC_DT -UMD_EDHEFLEX_WALL -UMD_INELASTIC -UMD_POLYDISP -UMD_CALC_DPP -UMD_HE_PARALL -DMD_BIG_DT -DMD_ASYM_ITENS -DEDHE_FLEX -DMD_PATCHY_HE -UMD_GLOBALNR -DMD_GLOBALNRD -DMD_GLOBALNRNL -DMD_GLOBALNRDNL -UMD_USE_CBLAS -UMD_USE_LAPACK -DMD_STOREMGL -UMD_HSVISCO -DMD_LOADMESH -DXDISP -UATPTENS -DMOLPTENS -UATPRESS' #getmacros()
    MACROS=getmacros()
    flags = [ '-x', 'c', '-Wall', '-Wextra', '-I', DIR_OF_THIS_SCRIPT+'/commSrc']
    #flags.append(MACROS)
    flags=flags+MACROS
    return {'flags': flags}
