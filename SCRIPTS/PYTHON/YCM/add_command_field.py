#!/usr/bin/env python3
import sys
fn='compile_commands.json'
with open(fn) as f:
    lines=f.readlines()
# definisco un tipo compatibile con quello che si trova nei json file
# ossia una lista di dizionari
#questo serve per evitare degli errori poichè la json1 e json1 vengono create a runtime
#con la exec
json=[{'a':['']}]
jout=[{'a':['']}]
# creo due stringhe
stri=''.join(lines).replace('\n','')
#aggiungo le definizioni delle due variabili 
stri = 'json =' + stri
#eseguo str1b e str2b come codice python
exec(stri)
#scorri tutti i dizionari nella lista e ne faccio un merge tramite la funzione update che viene messo in json1
for c in json:
    a=' '.join(c['arguments'])
    ad = {'command': a}
    #c.update(ad)   
    print (ad,a)
    jout.append(c)
#stampa il risulta il risuultante json file
#print (jout)
