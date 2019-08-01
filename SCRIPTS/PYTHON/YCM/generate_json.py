#!/usr/bin/env python3
import json, os, sys
fn='compile_commands.json'
os.system('make clean')
args=sys.argv
if len(args) > 2:
    print('You have to supply one argument which will be passed to make')
    quit()
if len(args)==1:
    targ=''
else:
    targ = args[1]
#esegue la make tramite intercept-build e passa come argomento di make 
#l'argomento fornito al presente script
os.system('intercept-build make '+targ)
with open(fn) as f:
    jf=json.load(f)
# jf è una lista di dizionari
# li scorro tutti e aggiungo ad ogni dizionario
# il campo 'command' che si ottiene facendo il join in unica
# stringa del campo arguments in cui ci sono tutti gli argomenti
# ovvero se ho il comando 'gcc -Wall file.c'
# 'arguments' : ['gcc','-Wall', 'file.c']
# mentre 'command' sarà
# 'command': 'gcc -Wall file.c'
for i in range(0,len(jf)):
    args=jf[i]['arguments']
    command=' '.join(args) 
    jf[i]['command'] = command
#stampa il risultante json file in maniera leggibile, i.e.indentando e andando a capo 
#prima lo metto in una stringa
jfs = json.dumps(jf, indent=4)
#poi stampo la stringa su disco sovrasrivendo il file json
with open(fn,'w+') as f:
    f.write(jfs)
