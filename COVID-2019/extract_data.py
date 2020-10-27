#!/usr/local/bin/python3
import sys,json,datetime
#defaults
provincia='Terni'
Regione='Umbria'
datatype='n'
args=sys.argv

#helper functions
def print_error():
    print('syntax: extract_data.py [--start|-s <giorno> ] [--exclude|-e <province o regioni da escludere separate da virgola>|-r|--regioni')
    print('|-p|--province|-n|--nazionale] [--relative|-rel] <campo (deceduti,terapia_intensiva,totale_casi,...)>')
    print('nota: per i possibili campi vedere i file json')

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

#args parsing
if len(args)==1:
    print_error()
    quit()
itargs=iter(args)
toexcl=[]
toincl=[]
firstday=0
printrel=False
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
    elif a == '-s' or  a == '--start':
        firstday = int(next(itargs))
    elif a == '-rel' or a== '--relative':
        printrel=True
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

#build filename with datasets
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

#read json file (it will be stored in a python dictionary)
fn = fn + '.json'
with open(fn) as f:
    # jf è una lista di dizionari
    jf=json.load(f)

#extract data according to selected field
valarr=[int]
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

#print out cumlative data
g = 0
gmax = len(valarr)
print('{',end='')
for v in valarr:
    if g < firstday:
        g+=1
        continue
    if printrel:
        va = valarr[g] - valarr[firstday]
    else:
        va = valarr[g]
    print ('{',g-firstday, ',', va,'}',end='')
    if g < gmax-1:
        print(',',end='')
    g=g+1
print('}')
