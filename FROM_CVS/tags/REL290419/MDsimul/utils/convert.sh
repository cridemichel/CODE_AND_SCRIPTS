#!/bin/sh
for FF in `ls Cnf*_R*` 
do
$HOME/MDsimul/bin/cri2fraMIX $FF ${FF}_FRA
done
