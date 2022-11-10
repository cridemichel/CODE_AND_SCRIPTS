#!/bin/bash 

function test_and_run {
    echo -n $1:
    ping -q -c 1 $1 > /dev/null
    echo $?
    ssh $1 "rm -f /local/studente/*.{c,c~,C,C~,out,exe,dat,py,py~,o,sh,sh~,java,java~,pdf,txt}"
    ssh $1 "find /local/studente -iname \*.c -exec rm -rf {} \;"
    ssh $1 "find /local/studente -iname \*.c# -exec rm -rf {} \;"
    ssh $1 "find /local/studente -iname \*.dat -exec rm -rf {} \;"
    ssh $1 "find /local/studente -iname \*~ -exec rm -rf {} \;"
    ssh $1 "find /local/studente -iname \*.py -exec rm -rf {} \;"
    ssh $1 "find /local/studente -iname \*.x -exec rm -rf {} \;"
    ssh $1 "find /local/studente -iname \*.png -exec rm -rf {} \;"
    ssh $1 "rm -rf /local/studente/LCNG* /local/studente/LCCDM*"
}  

for j in `seq 1 41` 
do 
    test_and_run 192.168.1.$j
done
