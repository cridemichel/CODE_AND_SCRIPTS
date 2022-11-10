#!/bin/bash

#mkdir $dir
cc=0
for i in `seq 1 1 41`
do
    ping -q -W 1 -c 1 192.168.1.$i > /dev/null
    if [ $? -eq 0 ]
    then
        rsync -av --delete 192.168.1.$i:'/local/studente/LCCDM*' . &> _aaa_
        if [ -f _aaa_ ]
        then
           FOL=$(cat _aaa_ | tail -n4 | grep LCCDM | awk -F / '{print $1}')
        else
           FOL=""
        fi
        if [ "$FOL" != "" ] 
	then
          echo "192.168.1.$i $FOL"
          cc=$[$cc+1]
        else
          echo "192.168.1.$i has no LCCDM folder or nothing has been updated"
        fi
        rm _aaa_ > /dev/null
    fi
done

copied=$(ls -1drt LCCDM* | wc -l)
echo "$copied cartelle copiate"
rm -f esercitazioni.zip &> /dev/null
zip -r esercitazioni.zip LCCDM* &> /dev/null
chown -R demichele.docenti .
mv esercitazioni.zip /home/docenti/demichele/
