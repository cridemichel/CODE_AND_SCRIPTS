#!/bin/bash
# N.B. ricordarsi di commentare la riga
# export CC=$(GCC) nel Makefile
# prima di lanciare questo script
# inoltre vanno installati scan-build e compdb con i seguenti comandi
# pip3 install scan-build compdb
if [ "$1" != "" ]
then
  TARGET="$1"
else
  TARGET=""
fi
COMPDB="1"
if [ "$2" == "0" ]
then
  COMPDB="0"
fi
make clean
intercept-build make $TARGET
if [ "$COMPDB" == "1" ]
then
compdb -p . list > _bbb_
mv compile_commands.json _aaa_
./merge_json.py _aaa_ _bbb_ > compile_commands.json
#mv -f _bbb_ compile_commands.json
mv compile_commands.json _aaa_
#riformatta il file in maniera più leggibile
#jq si installa con 'brew install jq'
jq . _aaa_ > compile_commands.json
rm -f _bbb_
fi
