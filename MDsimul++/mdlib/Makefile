#commSrc Makefile
OBJS= prova.o
CXX=g++
LDFLAGS=-g
CFLAGS=-g --std=c++17
common: $(OBJS)

prova: prova.o
	$(CXX) $(LDFLAGS) prova.o -o prova -lm	

prova.o : prova.c++ rpoly.H pvector.H simulation.H pmatrix.H
	$(CXX) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o prova.o prova.c++ 

clean:
	rm -f *.o

realclean: 
	rm -f *.o 
