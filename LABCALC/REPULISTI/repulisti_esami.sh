#!/bin/bash 

function test_and_run {
    echo -n $1:
    ping -q -c 1 $1 > /dev/null
    echo $?
    ssh $1 "rm -rf /local/studente/EXLC_*"
}  

for j in `seq 1 41` 
do 
    test_and_run 192.168.1.$j
done
