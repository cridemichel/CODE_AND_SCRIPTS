#!/usr/local/bin/python3
import os, sys
arg=sys.argv
if len(arg) < 2:
    print ('You must supply a pdf file')
    print ('extract_workers.py <pdf_file>')
    quit()
FPDF=arg[1]
#NOMI=arg[2]
INC=2
i=1
c=1
ii=i+2
FTXT='testo.txt'
os.system("pdftotext '" + FPDF + "' " + FTXT)
with open(FTXT, encoding='utf-8') as f:
    lines=f.readlines()
#with open(NOMI) as f:
#    nomi=f.readlines()
np=0 #numero pagine 
nw=0 #numero lavoratori
pagini=1
listpag=[]
cc=0
nomi=[]
fine=False
tempwork=False
nomco=''
codfisc=''
cf_changed=False
fine=False
found_nc = False
last_found_nc=1
for l in lines:
    #print ('linea=',l)
    if  l.find('Zucchetti spa, Autorizzazione') != -1:
        np = np+1
        # qui siamo alla fine della pagina, quindi alla fine della prima pagina np=1
        # alla fine della seconda np=2, ecc.
        if cf_changed==True:
            # se il codice fiscale è cambiato nella pagina attuale, si tratta di un nuovo 
            # lavoratore quindi va aggiunto il suo nome alla lista dei nomi (che si chiama 'nomi')
            # e vanno salvate le pagine del lavoratore precedente 
            #print('cf_changed pg. N. ', np)
            cf_changed=False
            nw = nw+1
            sublines = lines[cc:]
            # il primo lavoratore ha come pagine iniziale 1 e come pagina finale
            # l'ultima pagina dove è stato trovato un nome e cognome prima di tale cambiamento
            if nw > 1:
                # la pagina finale del lavoratore precedente è l'ultima dove 
                # è stato trovato un nome
                pagfin=last_found_nc
                listpag.append([pagini,pagfin])
                # la pagina successiva all'ultima dove è stato trovato un nome è la
                # la prima del lavoratore successivo
                pagini=last_found_nc+1
            nomi.append(nomco)
        if fine==True:
            # siamo alla prima pagina con la stringa "RIEPILOGO" quindi abbiamo finito
            # la pagina precedente è l'ultima dell'ultimo lavoratore
            pagfin = np - 1
            print('pag fin=', pagfin)
            listpag.append([pagini,pagfin])
            break
        if found_nc == True:
            # memorizzo la pagina attuale poiché abbiamo trovato un nome e cognome
            last_found_nc = np
            found_nc=False
    # la linea con il RIEPILOGO è prima di quella con 'Zucchetti spa,...'
    # quindi np è la pagina precedente ovvero l'ultimna dell'ultimo lavoratore
    if l.find('RIEPILOGO') != -1:
        fine=True
    delcc=-1
    # assumo che la stringa COGNOMEsEsNOME sia presente nell'ultima
    # pagina del cedolino di ogni lavoratore
    if l.find('COGNOMEsEsNOME')!=-1 and lines[cc+4].find('RIEPILOGO')==-1:
        #print('l1=', l)
        found_nc = True
        nomcogn_found=True
        delcc = 6
        nomco_old = nomco
        codfisc_old = codfisc
        # attualmente il nome e cognome sono a +6 dalla riga con COGNOMEsEsNOME
        # e due righe più giù c'è il codice fiscale
        nomco = lines[cc+delcc].strip('\n').replace(" ","_").replace("'","_")
        codfisc = lines[cc+delcc+2].strip('\n')
        #print(nomco, '<>', nomco_old)
    if delcc!=-1 and codfisc != codfisc_old:
        cf_changed=True
    cc=cc+1
#listpag.append([pagini,np])
print ('# pagine=', np+1)
print ('# lavoratori=', nw, ' ( #PDF:', len(listpag), ' #NOMI:', len(nomi), ')')
subdir='CEDOLINI'
if not os.path.exists(subdir):
    os.system('mkdir '+subdir)
#for n in nomi:
#    print('n=', n);
cc=0
print('nomi=',nomi)
for n in nomi:
    pag=listpag[cc]
    print(nomi[cc].strip('\n') + ' pagine da ' + str(pag[0]) + ' a ' + str(pag[1]))
    os.system('gs -dSAFER -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFirstPage='+ str(pag[0]) + ' -dLastPage='+ str(pag[1]) +' -sOUTPUTFILE=' + nomi[cc].strip('\n') + '.pdf ' + "'"+FPDF+"'")
    #gs -dSAFER -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFirstPage=$i -dLastPage=$ii -sOUTPUTFILE=$NOME.pdf $FN
    os.system('mv ' + nomi[cc].strip('\n') + '.pdf ' + subdir + '/')
    cc=cc+1
os.system('rm ' + FTXT)
