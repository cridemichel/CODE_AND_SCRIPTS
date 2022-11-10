#!/bin/bash
#LS="/_COSTANTINI|_GREALEY|_PISATURO|_ROSSI|_PERALTA|_SELLESE|_MAJERCZYK|_CAPOGRASSI/{print $0}"
#mkdir $dir
cd /home/docenti/demichele/CDM/LC/
cc=0
for i in `seq 1 1 41`
do
    ping -q -W 1 -c 1 192.168.1.$i > /dev/null
    if [ $? -eq 0 ]
    then
        rsync -av --delete 192.168.1.$i:'/local/studente/*_*' . &> _aaa_
        if [ -f _aaa_ ]
        then
           FOL=$(cat _aaa_ | tail -n4 | awk '/_RINALDI|_COSTANTINI|_GREALEY|_PISATURO|_ROSSI|_PERALTA|_SELLESE|_MAJERCZYK|_CAPOGRASSI/{print $0}' | awk -F / '{print $1}')
        else
           FOL=""
        fi
        if [ "$FOL" != "" ] 
	then
          echo "192.168.1.$i $FOL"
          cc=$[$cc+1]
        else
          echo "192.168.1.$i has no student folder or nothing has been updated"
        fi
        rm _aaa_ > /dev/null
    fi
done

copied=$(ls -1drt *_* | wc -l)
echo "$copied cartelle copiate"
rm -f esame_labcalc_maggio22.zip &> /dev/null
zip -r esame_labcalc_maggio22.zip [A-Z]*_[A-Z]* &> /dev/null
chown -R demichele.docenti .
mv esame_labcalc_maggio22.zip /home/docenti/demichele/
