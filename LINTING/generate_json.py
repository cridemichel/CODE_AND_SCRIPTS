#!/usr/bin/env python3
import json, os, sys

#  To manually create a json file for a specific file by Clang
#  Clang's -MJ option generates a compilation database entry per input (requires Clang >= 5.0).
#  Usage:
#  clang++ -MJ a.o.json -Wall -std=c++11 -o a.o -c a.cpp
#  clang++ -MJ b.o.json -Wall -std=c++11 -o b.o -c b.cpp
#  To merge the compilation database entries into a valid compilation database, it is possible to use (GNU) sed:
#  sed -e '1s/^/[\n/' -e '$s/,$/\n]/' *.o.json > compile_commands.json
#  Or, using any sed under Bash, Zsh or ksh:
#  sed -e '1s/^/[\'$'\n''/' -e '$s/,$/\'$'\n'']/' *.o.json > compile_commands.json
#  This sed invocation does the following:
#      insert the opening bracket: [
#      concatenate the entries
#      remove the trailing comma of the last entry (to be JSON compliant)
#      insert the closing bracket: ]

fn='compile_commands.json'
os.system('make clean')
use_compdb=False
copy_json=False
args=sys.argv
del(args[0])
itargs=iter(args)
targ=''
engine='bear --' # default engine
neng=1
# each element is an engine name, its abbreviated version and the full command to use from the shell
englist=[['compiledb','cdb','compiledb make '],['bear', 'be', 'bear -- make CXX=clang++ CC=clang '],['intercept-build', 'ib', 'intercept-build make CXX=clang++ CC=clang'],['clang', 'cl','make CXX="clang++ -MJ -" CC="clang -MJ -"']]
for a in itargs:
    if a=='--headers' or a == '-ih':
        use_compdb=True
    elif a=='--engine' or a=='-e':
        try:
            b=next(itargs)
        except StopIteration:
            print('[ERROR] no engine supplied, please provide one of the following ones:')
            print("compiledb, bear, intercept-build or clang")
            print("As of 11/11/2021 intercept-build produces an empty json file")
            quit()
        engfound=False
        neng=0
        for en in englist:
            if b == en[0] or b == en[1]:
                engine = en[2]
                engfound=True
                break
            neng = neng + 1
        if not engfound:
            print("[ERROR] Wrong engine name, it should be one of the following ones:")
            print("compiledb, bear, intercept-build or clang")
            print("As of 11/11/2021 intercept-build produces an empty json file")
            quit()
    elif a=='--copy' or a=='-c':
        copy_json=True
    elif a=='-h' or a=='--help':
        print("generate_json.py [-c|-ih] <arguments passed to make>")
        print("-c: copy the compile_commands.json file to parent directory (..)")
        print("-ih: include header in the compile_commands.json file")
        print("--engine|-e <engine>: specifies the engine to use. i.e. compiledb, bear, intercept-build or clang")
        print("As of 11/11/2021 intercept-build produces an empty json file.")
        print("Bear engine must be used in conjuction with mac osx g++/gcc (clang)")
        quit()
    else:
        targ=targ+' '+a
#esegue la make tramite intercept-build e passa come argomento di make 
#l'argomento fornito al presente script
print("Using engine " + englist[neng][0] + " ...")
if neng == 1: # neng=1 is the bear engine (see above)
    print("Bear engine must be used in conjuction with mac osx g++/gcc (clang),")
    print("otherwise an empty json file is obtained")
if neng == 3:
    os.system(engine + targ + ' > _genjson_.out')
    with open('_genjson_.out', encoding='utf-8') as f:
        lines=f.readlines()
    nl=['[']
    for l in lines:
        lt=l.strip('\n')
        larr=lt.split()
        # skip lines which are the ones produced by -MJ flag
        # these lines begin with { and "directory":
        if larr[0]=='{' and larr[1]=='"directory":':
            nl.append(lt)
    # remove the comma (,) from last line
    s=nl[-1][:-1]
    nl[-1]=s
    nl.append(']')
    with open('compile_commands.json','w+',encoding='utf-8') as f:
        for l in nl:
            f.write(l)
    os.system('rm _genjson_.out')
else:
    os.system(engine + targ)
if use_compdb:
    # add header files to compile_commands.json
    # compdb read a compile_commands.json file and add entries
    # for all the headers in the project recursively!
    os.system('compdb list > _compile_commands.json_')
    os.system('mv _compile_commands.json_ compile_commands.json')
    # compdb generates only "command" entries but not "arguments" entries 
    # in the following we add them for all files
    with open(fn, encoding='utf-8') as f:
        jf=json.load(f)
    for i in range(0,len(jf)):
        commands=jf[i]['command']
        jf[i]['arguments']=commands.split()
    jfs = json.dumps(jf, indent=4)
    #poi stampo la stringa su disco sovrasrivendo il file json
    with open(fn,'w+',encoding='utf-8') as f:
        f.write(jfs)
#
else:
    with open(fn,encoding='utf-8') as f:
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
    with open(fn,'w+',encoding='utf-8') as f:
        f.write(jfs)
if copy_json==True:
    os.system('cp compile_commands.json ..')
