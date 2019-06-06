#!/usr/bin/env python3
import sys
import os
def get_proc_cmdlines(p=''):
	cls=[]
	allpids=[]
	pids = [pid for pid in os.listdir('/proc') if pid.isdigit()]
	for pid in pids:
		try:
			with open(os.path.join('/proc', pid, 'cmdline'), 'rb') as f:
				l=f.readline()
				l=l.replace(b'\x00',b'\x20')
				l2=str(l.decode('utf-8'))
			if p in l2:
				cls.append(l2)			
				allpids.append(pid)	
		except IOError: # proc has already terminated
			continue
	return (cls,allpids)
#get X0
args=sys.argv
if len(args) > 1:
	x0 = args[1]
else:
	dir=os.getcwd()
	a=dir.split('/')
	b=a[len(a)-1].split('_')
	if len(b) < 2:
		print('[ERROR] you did not supply an X0 and you are in the right place')
		quit()
	x0=b[1]
#print(x0)
c=0
#return all pid obtained from commandlines containing string p
l2,pids=get_proc_cmdlines('veff_IR_')
#print("pids=",pids)
allIR=[]
for e in l2:
	l=e.split(' ')
	if len(l) == 8:
		if l[3] == x0:
			#print("process",pids[c],"=",l)
			c += 1
			allIR.append(l[5])
#
#print ('running=', c)
#print ('allIR=', allIR)
#search for missing jobs and check if they have already finished
perc='/home/demichel/HARDELL/Veff/X0_'+str(x0)+'/'
ok=True
nd=0
if not os.path.exists(perc):
	print('dir \''+perc+'\' does not exist...\n')
	quit()
for i in range(0,300) :
	if str(i) not in allIR:
		#print ('job ', i, ' is missing')
		dir=perc+'IR_'+str(i)
		try:
			with open(dir+'/veff_vs_tt.dat') as f:
				lines=f.readlines()
				lastline=lines[-1]
				l=lastline.strip('\n').split(' ')
				if l[0] != '99900000000':
					print('job IR='+str(i)+' is not running and it has not finished yet!')
					print('I am restarting it')
					exec=dir+'/veff_IR_'+str(i)
					print ('exec is: '+exec)
					ok=False
				else:
					nd+=1
		except:
			print ('file '+dir+'/veff_vs_tt.dat'+ 'does not exist but I continue...\n')
		continue
if not ok:
	print('Some jobs were dead and I had to restart them!\n')
else:
	if c == 0:
		print('All done here!')
	else:
		print('Here (X0='+ x0 + ') there are '+str(c)+' jobs runnings', end='')
		print(' and '+ str(nd) + ' regularly finished')
		print('Total job is', str(nd+c))
