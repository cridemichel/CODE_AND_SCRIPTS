#!/usr/local/bin/python3
import sys,json,datetime
#from datetime import datetime
fn='compile_commands.json'
provincia='Terni'
Regione='Umbria'
args=sys.argv
def print_error():
    print('syntax: extract_data.py [--exclude|-e <province o regioni da escludere separate da virgola>|-r|--regioni')
    print('|-p|--province|-n|--nazionale] <campo (deceduti,terapia_intensiva,totale_casi,...)>')
    print('r=regiorni, p=province, n=nazionale')

def parse_ranges(arg):
    sarg=[]
    if arg.find(',')!=-1:
        sarg=arg.split(',')
    else:
        sarg.append(arg)
    return sarg

def data2obj(d):
    subs = d.split('T')[0]
    return datetime.datetime.strptime(subs,"%Y-%m-%d")
def inlist(n,l):
    if n in l:
        return True
    else:
        return False
if len(args)==1:
    print_error()
    quit()
itargs=iter(args)
toexcl=[]
toincl=[]
for a in itargs:
    if a == '--help' or a  == '-h':
        print_error()
        quit()
    elif a == '-i' or a == '--include':
        toincl = parse_ranges(next(itargs))
    elif a == '-e' or a == '--exclude':
        toexcl = parse_ranges(next(itargs))
    elif a == '-p' or a == '--province':
        datatype='p'
    elif a == '-r' or a == '--regioni':
        datatype='r'
    elif a == '-n' or a == '--nazionale':
        datatype = 'n'
    else:
        datafield = a
if len(toexcl) > 0 and len(toincl) > 0:
    print('You can use either exclude or include option')
    quit()
if len(toexcl) > 0:
    seltype='e'
elif len(toincl) > 0:
    seltype='i'
else:
    seltype='a' #all
fn = 'dpc-covid19-ita-'
if datatype == 'p':
    fn = fn + 'province'
elif datatype == 'r':
    fn = fn + 'regioni'
elif datatype == 'n':
    fn = fn + 'andamento-nazionale'
else:
    print_error()
    quit()
fn = fn + '.json'
with open(fn) as f:
    jf=json.load(f)
# jf è una lista di dizionari
valarr=[int]
#for i in range(0,100):
#    valarr.append(0)
valarr=[]
primo=True
for i in range(0,len(jf)):
    if datatype == 'p':
        provincia= jf[i]['denominazione_provincia']
        if seltype=='e':
            if inlist(provincia,toexcl):
                continue
        elif seltype=='i':
            if not inlist(provincia,toincl):
                continue
    elif datatype == 'r':
        regione = jf[i]['denominazione_regione']
        if seltype=='e':
            if inlist(regione,toexcl):
                continue
        elif seltype=='i':
            if not inlist(regione,toincl):
                continue
    if primo==True:
        dataini=data2obj(jf[i]['data'])
        primo=False
    datarel = data2obj(jf[i]['data'])-dataini
    valore=jf[i][datafield]
    giorno = datarel.days
    #giorno = int(numpy.rint((datarel.seconds)/(60*24)))
    #print('seconds=', datarel.seconds, ' giorni=', datarel.days)
    if giorno > len(valarr)-1:
        valarr.append(valore)
    else:
        valarr[giorno] += valore
g = 0
gmax = len(valarr)
print('{',end='')
for v in valarr:
    print ('{',g, ',', valarr[g],'}',end='')
    if g < gmax-1:
        print(',',end='')
    g=g+1
print('}')
