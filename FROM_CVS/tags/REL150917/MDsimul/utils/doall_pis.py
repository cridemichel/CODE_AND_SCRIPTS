#!/usr/bin/python2.2
#in lista ogni riga e' nel formato seguente:
#<directory> <temperatura>
#con un solo spazio che separi i due campi
import os
import sys
import string
fname = sys.argv[1]
os.system("rm pis-vs-T.dat")
F=file(fname,"r")
lines = F.readlines()
for ll in lines:
	print ll
	A = ll[0:(len(ll)-1)].split(" ")
	os.system("echo `pwd`/ > lista")
	os.system("ls " + A[0] + "Qnc* >> lista")
	os.system("/home/demichel/MDsimul/bin/calcpress lista | awk \'{ print $1,$2}\'> pis.tmp ")
	os.system("/home/demichel/MDsimul/utils/mean.sh pis.tmp | awk \'{print $2," + A[1] + "}\' >> pis-vs-T.dat ")
	os.system("rm pis.tmp");
