#!/bin/bash
for f in `cat lista`
do
   PVAL="$f"
   cat submit_$f | gawk -v pval=$PVAL -F '=' '{if ($1=="_PVAL_") printf("PVAL=\"%s\"\n", pval); else print($0);}'> _aaa_
   cat _aaa_     | gawk -v jname="univscal_$PVAL" -F '=' '{if ($2=="_JNAME_") printf("%s=%s\n", $1, jname); else print($0);}'> _bbb_
   rm _aaa_
   mv _bbb_ submit_$f   
done
