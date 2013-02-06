if [ \( "$1" == "" \) -o \( "$2" == "" \) ]
then
echo "ERROR: kill_sus_range.sh <min N> <max N>"
exit
fi
ps ax |grep HC| awk -F _ -v MAXN=$2 -v MINN=$1 '{if ($3 < MAXN && $3 > MINN) print $0}' | awk '{print("kill", $1)}' > _aaa_
sh _aaa_
rm _aaa_
