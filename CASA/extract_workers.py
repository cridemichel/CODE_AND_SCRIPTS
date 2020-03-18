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
os.system('pdftotext ' + FPDF + ' ' + FTXT)
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
for l in lines:
    #print ('linea=',l)
    if  l.find('VIDIMAZIONE') != -1:
        np = np+1
    if  l.find('AGGIUNTIVE') != -1:
        nw = nw+1
        pagfin = np-1
        if nw > 1:
            listpag.append([pagini,pagfin])
        pagini = np
    if l.find('COGNOME') != -1:
        sublines = lines[cc:]
        ccc=0
        for sl in sublines:
            #print('subline=', sl)
            if sl.find('SALUS')!=-1:
                nomco=sublines[ccc+1]
                nomus=nomco.replace(' ', '_')
                nomi.append(nomus)
                break
            ccc = ccc+1
        #print('nome cognome=', nomco)
    cc=cc+1
listpag.append([pagini,np-3])
print ('# pagine=', np)
print ('# lavoratori=', nw)
for n in nomi:
    print('n=', n);
cc=0
for n in nomi:
    pag=listpag[cc]
    print('pag=', pag[0], pag[1])
    os.system('gs -dSAFER -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFirstPage='+ str(pag[0]) + ' -dLastPage='+ str(pag[1]) +' -sOUTPUTFILE=' + nomi[cc].strip('\n') + '.pdf ' + FPDF)
    #gs -dSAFER -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFirstPage=$i -dLastPage=$ii -sOUTPUTFILE=$NOME.pdf $FN
    os.system('mv ' + nomi[cc].strip('\n') + '.pdf SUBCEDOL/')
    cc=cc+1
os.system('rm ' + FTXT)
