with open('hosts_1',encoding='utf8') as f:
    al=f.readlines()
fl = []
for l in al:
    sl=l.strip('\n').split()
    if sl=='\n':
        continue
    if len(sl) < 2:
        continue
    if not (sl[1] in fl):
        fl.append(sl)
for l in fl:
    print(l[0],' ', l[1])
