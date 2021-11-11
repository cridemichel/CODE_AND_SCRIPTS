#!/usr/bin/env python3
import json, os, sys
fn='compile_commands.json'
os.system('make clean')
use_compdb=False
copy_json=False
args=sys.argv
del(args[0])
itargs=iter(args)
targ=''
engine='bear --' # default engine
# each element is an engine name, its abbreviated version and the full command to use from the shell
englist=[['compiledb','cdb','compiledb'],['bear', 'b', 'bear --'],['intercept-build', 'ib', 'intercept-build']]
for a in itargs:
    if a=='--headers' or a == '-ih':
        use_compdb=True
    elif a=='--engine' or a=='-e':
        b=next(itargs)
        engfound=False
        for en in englist:
            if b == en[0] or b == en[1]:
                engine = en[2]
                engfound=True
        if not engfound:
            print("[ERROR] Wrong engine name, it should be one of the following ones:")
            print("compiledb, bear or intercept-build")
            print("As of 11/11/2021 intercept-build produces an empty json file")
            quit()
    elif a=='--copy' or a=='-c':
        copy_json=True
    elif a=='-h' or a=='--help':
        print("generate_json.py [-c|-ih] <arguments passed to make>")
        print("-c: copy the compile_commands.json file to parent directory (..)")
        print("-ih: include header in the compile_commands.json file")
        print("--engine|-e <engine>: specifies the engine to use. i.e. compiledb, bear or intercept-build")
        print("As of 11/11/2021 intercept-build produces an empty json file")
        quit()
    else:
        targ=targ+' '+a
#esegue la make tramite intercept-build e passa come argomento di make 
#l'argomento fornito al presente script
os.system(engine + ' make ' + targ)
if use_compdb:
    # add header files to compile_commands.json
    # compdb read a compile_commands.json file and add entries
    # for all the headers in the project recursively!
    os.system('compdb list > _compile_commands.json_')
    os.system('mv _compile_commands.json_ compile_commands.json')
    # compdb generates only "command" entries but not "arguments" entries 
    # in the following we add them for all files
    with open(fn) as f:
        jf=json.load(f)
    for i in range(0,len(jf)):
        commands=jf[i]['command']
        jf[i]['arguments']=commands.split()
    jfs = json.dumps(jf, indent=4)
    #poi stampo la stringa su disco sovrasrivendo il file json
    with open(fn,'w+') as f:
        f.write(jfs)
#
else:
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
if copy_json==True:
    os.system('cp compile_commands.json ..')
