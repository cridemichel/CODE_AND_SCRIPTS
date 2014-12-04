#!/bin/sh
cat $1 |  awk -v NUM=$2 '{\
 if ($1 == "#define" && $2 == "MD_SIMDAT")\
   {\
     if (NUM==0)\
       print ($1 " " $2 " " $3 " " "\"/simdat\"" );\
     else\
       print ($1 " " $2 " " $3 " " "\"/simdat" NUM "\"" );\
   }\
 else print $0   }' > _AWKOUT_
 DIFF=`diff _AWKOUT_ $1`
 if [ "$DIFF" != "" ] 
 then
  mv -f _AWKOUT_ $1
 else
  rm -f _AWKOUT_
 fi
