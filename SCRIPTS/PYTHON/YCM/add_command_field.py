#!/usr/bin/env python3
import json
fn='compile_commands.json'
with open(fn) as f:
    jf=json.load(f)
#scorri tutti i dizionari nella lista e ne faccio un merge tramite la funzione update che viene messo in json1
for i in range(0,len(jf)):
    args=jf[i]['arguments']
    command=' '.join(args) 
    jf[i]['command'] = command
#stampa il risultante json file in maniera leggibile, i.e.indentando e andando a capo 
jfs = json.dumps(jf, indent=4)
with open(fn,'w+') as f:
    f.write(jfs)
