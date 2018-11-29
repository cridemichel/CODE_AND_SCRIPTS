#commSrc Makefile
OBJS= prova.o
CXX=g++
LDFLAGS=-g
OPTS=-ffast-math -O3 -march=native -mtune=native -mavx2 #-march=native sul laptop migliora di un 30% le prestazioni! 
CFLAGS=-g --std=c++17 # $OPTS
common: $(OBJS)

prova: prova.o
	$(CXX) $(LDFLAGS) prova.o -o prova -lm	

prova.o : prova.c++ rpoly.H pvector.H simulation.H pmatrix.H
	$(CXX) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o prova.o prova.c++ 

clean:
	rm -f *.o

realclean: 
	rm -f *.o 
