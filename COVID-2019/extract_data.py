#!/usr/local/bin/python3
import sys,json,datetime
from datetime import datetime
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
    return datetime.strptime(d,"%Y-%m-%dT%H:%M:%S")

if len(args)==1:
    print_error()
    quit()
itargs=iter(args)
toexcl=[]
for a in itargs:
    if a == '--help' or a  == '-h':
        print_error()
        quit()
    elif a == '-kill' or a == '-k':
        toexcl=parse_ranges(next(itargs))
    elif a == '-p' or a == '--province':
        datatype='p'
    elif a == '-r' or a == '--regioni':
        datatype='r'
    elif a == 'n' or a == '--nazionale':
        datatype = 'n'
    else:
        datafield = a
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
for i in range(0,len(jf)):
    if i==0:
        dataini = data2obj(jf[i]['data'])
    datarel = data2obj(jf[i]['data'])-dataini
    valore=jf[i][datafield]
    giorno = int(datarel.days)
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
