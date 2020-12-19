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
tempwork=False
lastcf=-1
lastcn=-1
for l in lines:
    #print ('linea=',l)
    if  l.find('VIDIMAZIONE') != -1:
        np = np+1
    if  len(l) > 11 and l[0:6].isalpha() and l[6:11].isdigit():
        nw = nw+1
        pagini= np
        pagfin = np
        listpag.append([pagini,pagfin])
        nomco=l[11:]
        nomus=nomco.replace(' ', '_')
        nomi.append(nomus)
        lastcf=-1
    cc=cc+1
listpag.append([pagini,np])
print ('# pagine=', np)
print ('# lavoratori=', nw)
subdir='CEDOLINI'
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
