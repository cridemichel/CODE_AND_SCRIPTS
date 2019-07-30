#!/usr/bin/env python3
import os, time, sys
curt=time.time()
args=sys.argv
undodir='/Users/demichel/.vim/undo/'
if len(args) > 1:
    number_of_days=int(args[1])
else:
    number_of_days=180
def days_to_sec(d):
    return d*3600*24
for f in os.listdir(undodir):
    mtime=os.path.getmtime(undodir+'/'+f)
    if curt-mtime > days_to_sec(number_of_days):
        os.remove(undodir+'/'+f)
