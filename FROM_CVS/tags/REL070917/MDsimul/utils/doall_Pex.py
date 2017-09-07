#!/usr/bin/python2.2
#in lista ogni riga e' nel formato seguente:
#<directory> <temperatura>
#con un solo spazio che separi i due campi
import os
import sys
import string
fname = sys.argv[1]
F=file(fname,"r")
os.system("rm PexvsV.dat")
lines = F.readlines()
for ll in lines:
	print ll
	A = ll[0:(len(ll)-1)].split(" ")
	os.chdir(A[0])
	os.system("$HOME/MDsimul/bin/md2ascii -p 15 -o Pex-0.dat Pex-0");
	os.system("../mean.sh Pex-0.dat | awk -v VOL=" + A[1] + " \'{print VOL,$2}\' >> ../PexvsV.dat")
	os.system("cat PexvsV.dat | sort -n")
	os.chdir("..")
