#!/usr/local/bin/python3
import os, sys
arg=sys.argv
if len(arg) < 2:
    print ('You must supply a pdf file')
    print ('extract_workers_cud.py <pdf_file>')
    quit()
FPDF=arg[1]
#NOMI=arg[2]
INC=2
i=1
c=1
ii=i+2
FTXT='testo.txt'
os.system("pdftotext -layout '" + FPDF + "' " + FTXT)
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
tempwork=False
for l in lines:
    if  l.find('Pag.')!=-1:
        np = np + 1
        if l.find('1/4')!= -1:
            pagini = np
        if l.find('3/4')!= -1:
            pagfin = np
        if nw > 1:
            listpag.append([pagini,pagfin])
    if l.find('PENSIONATO') != -1:
        nw = nw + 1
        sublines = l.split()
        nomi.append(sublines[3]+'_'+sublines[4])
        #print('nome cognome=', nomco)
    cc=cc+1
print ('# pagine=', np)
print ('# lavoratori=', nw)
subdir='CUD'
if not os.path.exists(subdir):
    os.system('mkdir '+subdir)
#for n in nomi:
#    print('n=', n);
cc=0
for n in nomi:
    pag=listpag[cc]
    print(nomi[cc].strip('\n') + ' pagine da ' + str(pag[0]) + ' a ' + str(pag[1]))
    os.system('gs -dSAFER -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFirstPage='+ str(pag[0]) + ' -dLastPage='+ str(pag[1]) +' -sOUTPUTFILE=' + nomi[cc].strip('\n') + '.pdf ' + "'"+FPDF+"'")
    #gs -dSAFER -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFirstPage=$i -dLastPage=$ii -sOUTPUTFILE=$NOME.pdf $FN
    os.system('mv ' + nomi[cc].strip('\n') + '.pdf ' + subdir + '/')
    cc=cc+1
os.system('rm ' + FTXT)
