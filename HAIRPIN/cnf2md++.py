#!/usr/bin/python3
import sys,os
if len(sys.argv) < 2:
    print('You have to supply a filename')
    quit
fn=sys.argv[1]
with open(fn) as f:
    lines=f.readlines()
nat=0
nl=0
for l in lines:
    if l=='@@@\n':
        nat+=1;
        nl+=1
        continue
    if nat < 1:
        parname=l.split(':')[0].strip('\n')
        parval=l.split(':')[1].strip('\n')
#       print('parname=', parname, ' parval=', parval)
        if parname=='parnum':
            N=int(parval)
    elif nat==2:
        nbeg=nl
        break
    nl+=1
ll=lines[-1]
Lx=ll.strip('\n').split(' ')[0]
Ly=ll.strip('\n').split(' ')[1]
Lz=ll.strip('\n').split(' ')[2]
nend=len(lines)-2
print(N, ' ', Lx, ' ', Ly, ' ', Lz, ' 1 16.24 8 1 1 0')
print('0.5 1')
#print('nbeg=',nbeg, ' nend=', nend)
for l in lines[nbeg:nend+1]:
    nl=l.strip('\n').split(' ')[0:-1]
    print(' '.join(nl))
