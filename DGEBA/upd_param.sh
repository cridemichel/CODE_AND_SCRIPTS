#!/bin/bash
if [ "$1" == "" ] 
then
 echo "You must supply parameter name and value"
 exit 1
fi
if [ "$2" == "" ] 
then
 echo "You must supply parameter name and value"
 exit 1
fi
cat ellipsoid_dGEBA.par | awk -v PN=$1 -v PV=$2 -F ':' '{if ($1 == PN) printf("%s: %s\n", PN, PV); else print $0}' > _aaa_
mv _aaa_ ellipsoid_dGEBA.par
