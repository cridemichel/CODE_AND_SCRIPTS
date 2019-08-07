#!/usr/bin/env /usr/local/bin/python3
# USAGE:
# cleanupold.py <dir> <days> <regexp>
# delete all file in the directory <dir> older than <days> days
# if <regexp> is provided it deletes only files that matches the <regexp> provided
# put this into crontab for periodic cleanup
import os, time, sys, re
curt=time.time()
args=sys.argv
cleandir='/Users/demichel/.vim/undo/' #dir to cleanup
number_of_days=180
regexp=False # i.e. by default delete all files older than 'number_of_days' days
if len(args) > 1:
    cleandir= args[1]
if len(args) > 2:
    number_of_days=int(args[2])
if len(args) > 3:
    regexp=True
    pattern = re.compile(args[3])
else:
    regexp=False
def days_to_sec(d):
    return float(d*3600*24)
for f in os.listdir(cleandir):
    if regexp == True:
        if pattern.search(f):
            delfile=True
            #print ('cancello il file', f)
        else:
            delfile=False
    else:
        delfile=True
    mtime=os.path.getmtime(cleandir+'/'+f)
    #if delfile:
    #    print('file=', f, ' sec=', curt-mtime, ' days2sec=', days_to_sec(number_of_days), ' delFile=', delfile)
    if (delfile==True) and (curt-mtime > days_to_sec(number_of_days)):
        os.remove(cleandir+'/'+f)
