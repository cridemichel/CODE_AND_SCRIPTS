#!/usr/bin/python2.2
#in lista ogni riga e' nel formato seguente:
#<directory> <temperatura>
#con un solo spazio che separi i due campi
import os
import sys
import string
fname = sys.argv[1]
os.system("rm eis-e-vs-T.dat")
F=file(fname,"r")
lines = F.readlines()
for ll in lines:
	print ll
	A = ll[0:(len(ll)-1)].split(" ")
	os.system("cat " + A[0] + "Ene* > eneall.dat")
	os.system("mean.sh eneall.dat | awk \'{ print \"" + A[1] + "\",$2,$3}\' >> eis-e-vs-T.dat")
