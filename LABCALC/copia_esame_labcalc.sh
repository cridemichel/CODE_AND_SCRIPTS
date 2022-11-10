#!/bin/bash
#LS="/_COSTANTINI|_GREALEY|_PISATURO|_ROSSI|_PERALTA|_SELLESE|_MAJERCZYK|_CAPOGRASSI/{print $0}"
#mkdir $dir
cd /home/docenti/demichele/CDM/LC/
cc=0
for i in `seq 1 1 41`
do
    ping -q -W 1 -c 1 192.168.1.$i &> /dev/null
    if [ $? -eq 0 ]
    then
        rsync -av 192.168.1.$i:'/local/studente/EXLC_*' . &> _aaa_ 
        FOL=$(cat _aaa_ | grep EXLCL22)
        if [ "$FOL" != "" ] 
	then
          echo "192.168.1.$i $FOL"
          cc=$[$cc+1]
        else
          echo "192.168.1.$i has no student folder or nothing has been updated"
        fi
        rm _aaa_ &> /dev/null
    fi
done

copied=$(ls -1drt EXLCL22* | wc -l)
echo "$copied cartelle copiate"
rm -f esame_labcalc_novembre22.zip &> /dev/null
zip -r esame_labcalc_novembre22.zip EXLC_* &> /dev/null
cp esame_labcalc_novembre22.zip esame_labcalc$(date +'_%F_%T').zip
chown -R demichele.docenti .
cp esame_labcalc_novembre22.zip /root/
mv esame_labcalc_novembre22.zip /home/docenti/demichele/
