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
with open(FTXT) as f:
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
for l in lines:
    #print ('linea=',l)
    if  l.find('Zucchetti spa, Autorizzazione') != -1:
        np = np+1
    # la linea con il RIEPILOGO è prima di quella con 'Zucchetti spa,...'
    # quindi np è la pagina precedente ovvero l'ultimna dell'ultimo lavoratore
    if l.find('RIEPILOGO') != -1:
        pagfin = np
        print('pag fin=', pagfin)
        listpag.append([pagini,pagfin])
        break
#   if  l.find('INFORMAZIONI AGGIUNTIVE') != -1:
#        nw = nw+1
#        pagfin = np-1
#        if nw > 1:
#            listpag.append([pagini,pagfin])
#        pagini = np
    delcc=-1
    if l.find('COGNOMEsEsNOME')!=-1 and lines[cc+4].find('RIEPILOGO')==-1:
        #print('l1=', l)
        delcc = 6
        nomco_old = nomco
        codfisc_old = codfisc
        nomco = lines[cc+delcc].strip('\n').replace(" ","_").replace("'","_")
        codfisc = lines[cc+delcc+2].strip('\n').replace(" ","_").replace("'","_")
        #print(nomco, '<>', nomco_old)
    if delcc!=-1 and codfisc != codfisc_old:
        nw = nw+1
        pagfin=np-1
        sublines = lines[cc:]
        if nw > 1:
            listpag.append([pagini,pagfin])
        #ccc=0
        pagini=np
        #nomco=lines[cc+1].replace("'","_")
        #nomus=nomco.replace(' ', '_')
        nomi.append(nomco)
        #  for sl in sublines:
        #      #print('subline=', sl)
        #      # dopo aver trovato la string COGNOME
        #      # cerca l'occorenza della string SALUS
        #      # poiché nella linea successiva ci saranno
        #      # nome e cognome
        #      if sl.find('CodicesFiscale')!=-1:
        #          nomco=sublines[ccc+2].replace("'","_")
        #          nomus=nomco.replace(' ', '_')
        #          nomi.append(nomus)
        #          break
        #      ccc = ccc+1
        #  #print('nome cognome=', nomco)
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
