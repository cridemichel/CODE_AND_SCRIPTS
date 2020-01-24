#!/bin/bash
if [ "$1" == "" ]
then
  echo "Fornire il nome del file da controllare"
  echo "code_similarity_check.sh <nome_file>"
  exit
fi
jplag -l c/c++ -c $@
