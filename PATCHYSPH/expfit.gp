f(x)=A*x+B
fit f(x) 'persist.dat' using 1:(log($2)) via A,B
print "persistence length is ", -1/A
