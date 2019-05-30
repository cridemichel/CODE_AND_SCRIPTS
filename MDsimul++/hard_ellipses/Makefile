CXX=g++-9
CC=gcc-9
CXXFLAGS=-O3 -Wall -std=c++17#-g 
CFLAGS=-O3 #-Wall
MDLIB=../mdlib
HEADERS=$(MDLIB)/pvector.H $(MDLIB)/pmatrix.H $(MDLIB)/rpoly.H $(MDLIB)/boxes.H ./mcsim_hardell.H 
LDFLAGS=-lm -llapack -lblas #-g

all: mcsim veff create_triang_lattice expand

create_triang_lattice: create_triang_lattice.cpp $(HEADERS) $(MDLIB)/linked_cell_lists.H
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o create_triang_lattice create_triang_lattice.cpp 
	
expand: expand.cpp $(HEADERS) $(MDLIB)/linked_cell_lists.H
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o expand expand.cpp 
	
mcsim: mcsim_hardell.cpp $(HEADERS) $(MDLIB)/linked_cell_lists.H
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o mcsim_hardell mcsim_hardell.cpp 

ellipses_c: ellipses.c
	$(CC) $(CFLAGS) -c ellipses.c 

ellipses_cpp: mcsim_hardell.cpp $(HEADERS) $(MDLIB)/linked_cell_lists.H 
	$(CXX) $(CXXFLAGS) -DUSE_C_ROUTINES -c mcsim_hardell.cpp 

mcsimc: ellipses_c ellipses_cpp 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o mcsim_hardell mcsim_hardell.o ellipses.o

veff: calc_veff_hardell.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o calc_veff_hardell calc_veff_hardell.cpp 
#gcc-9 -c ellipses.c -O3
#g++-9 -lm -llapack -lblas prova_ellipses.o  ellipses.o -o testhe
clean:
	rm -fr calc_veff_hardell mcsim_hardell *.o
