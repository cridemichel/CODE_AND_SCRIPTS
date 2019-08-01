#!/usr/bin/env python3
import sys, os
args=sys.argv
if len(args) < 2:
    print ('you have two supply the names of two json files')
#leggo i due json file
with open(args[1]) as f1:
    lines1=f1.readlines()
with open(args[2]) as f2:
    lines2=f2.readlines()
# definisco un tipo compatibile con quello che si trova nei json file
# ossia una lista di dizionari
#questo serve per evitare degli errori poichè la json1 e json1 vengono create a runtime
#con la exec
json1=[{'a':['']}]
json2=[{'a':['']}]
# creo due stringhe
str1=''.join(lines1).replace('\n','')
str2=''.join(lines2).replace('\n','')
#aggiungo le definizioni delle due variabili 
str1b = 'json1 =' + str1
str2b = 'json2 =' + str2
#eseguo str1b e str2b come codice python
exec(str1b)
exec(str2b)
#scorri tutti i dizionari nella lista e ne faccio un merge tramite la funzione update che viene messo in json1
for i in range(0,len(json1)):
    json1[i].update(json2[i]) 
#stampa il risulta il risuultante json file
print(str(json1).replace("'", '"'))
