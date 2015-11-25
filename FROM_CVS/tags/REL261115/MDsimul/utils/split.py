#!/usr/bin/python
import sys
import os
import string
blen=100
outhdr=0
curperc=os.getcwd()
fname="lista"
eseguibile="basins"
#fname=file da splittare
if len(sys.argv) > 1:
	fname=sys.argv[1]
#lunghezza del singolo blocco (se < 0 => mette il percorso all'inizio del file 	
if len(sys.argv) > 2:
	blen=int(sys.argv[2])
if len(sys.argv) > 3:
	templatefile=sys.argv[3]
else:
	templatefile=""
blocco=1
def outscript(templatefile, blocco):
	if templatefile != "":
		S=file(templatefile, "r")
		templ = S.readlines()
		for ii in range(len(templ)):
			templ[ii]=str.replace(templ[ii], "@PWD@", curperc)
			templ[ii]=str.replace(templ[ii], "@EXEC@", eseguibile)
			templ[ii]=str.replace(templ[ii], "@LISTA@", sublista)
			templ[ii]=str.replace(templ[ii], "@HOME@", os.getenv("HOME"))
			templ[ii]=str.replace(templ[ii], "@TEMPLFILE@", templatefile+"."+str(blocco))
		S.close()
		S=file(templatefile+(".%d" % blocco),"w")
		S.writelines(templ)
		S.close()
		os.system("chmod a+x "+templatefile+"."+str(blocco))	
F=file(fname,"r")
sublista=fname+"."+str(blocco)
O=file(sublista,"w")
if blen < 0:
	blen = -blen
	outhdr=1;
	O.write(curperc+"/\n")
outscript(templatefile, blocco)
lista=F.readlines()
ll=len(lista)
i=0
while i < ll:
	O.write(lista[i])	
	i=i+1
	if i == (blocco*blen) and i < ll:
		blocco=blocco+1
		O.close()
		sublista = fname+"."+str(blocco)
		O = file(sublista,"w")
		if outhdr:
			O.write(curperc+"/\n")
		outscript(templatefile, blocco)	
