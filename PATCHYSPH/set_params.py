#!/usr/bin/python2.2
#in lista ogni riga e' nel formato seguente:
#<directory> <temperatura>
#con un solo spazio che separi i due campi
import os
import sys
import string
cc=2
while cc <= len(sys.argv)-1:
	#print(sys.argv[cc],sys.argv[cc+1])
	cmd="../set_one_param.sh "+sys.argv[1]+" "+sys.argv[cc]+" " +sys.argv[cc+1]
	#	print("cmd:"+cmd)
	os.system(cmd)
	cc=cc + 2
