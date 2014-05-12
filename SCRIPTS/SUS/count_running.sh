ps ax | grep SWHE_N | grep mosrun > _aaa_
wc -l _aaa_ | awk '{print $1}'
rm _aaa_
