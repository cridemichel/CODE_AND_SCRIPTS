#!/bin/bash
# N.B. ricordarsi di commentare la riga
# export CC=$(GCC) nel Makefile
# prima di lanciare questo script
# inoltre vanno installati scan-build e compdb con i seguenti comandi
# pip3 install scan-build compdb
if [ "$1" != "" ]
then
  TARGET="$1"
fi
make clean
intercept-build make $TARGET
compdb -p . list > _aaa_
mv -f _aaa_ compile_commands.json
