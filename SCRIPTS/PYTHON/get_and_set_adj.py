#!/usr/bin/python3
# DESCRIPTION:
# read values of deltra, vmax and delrot and set them in the xasc file supplied as argument
import sys
#curd=os.getcwd()
#d=os.path.split(curd.strip('\n'))[1]
args=sys.argv
lines2new=[]
with open('restart-0') as f:
    lines=f.readlines()
#print(lines)
deltra=lines[3].strip('\n').split()[2]
vmax=lines[3].strip('\n').split()[3]
delrot=lines[8].strip('\n')
print('Setting following values in ', args[1],': deltra=', deltra, ' vmax=', vmax, ' delrot=', delrot)
with open(args[1]) as f:
    lines2=f.readlines()
    for l2 in lines2:
        #print ('l2=', l2)
        ls=l2.strip('\n').split(':')
        p=ls[0].replace(' ', '')
        if p=='deltra':
            ln = 'deltra:'+deltra+'\n'
            lines2new.append(ln)
        #ln = ln.replace('_A_',sa)
        elif p=='vmax':
            ln = 'vmax:'+vmax+'\n'
            lines2new.append(ln)
        elif p=='delrot':
            ln = 'T:'+delrot+'\n'
            lines2new.append(ln)
        else:
            lines2new.append(l2)
with open(args[1],"w") as f:
    for l in lines2new:
        #print ('l=',l)
        f.write(l)
