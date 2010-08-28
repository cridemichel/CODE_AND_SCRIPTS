UGFACT="0.6"
UG=`echo "${UGFACT}*0.129^2" | bc -l`
#per awk log(x) = logaritmo naturale di x
cat $1 | LC_NUMERIC=C awk -v ug2=${UG} '{print (log(1.0/$1)/log(10.0)+log(ug2)/log(10.0),log($2)/log(10))}'
