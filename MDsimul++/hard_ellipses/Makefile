CXX=g++-9
CXXFLAGS=-O3 -Wall -std=c++17#-g 
MDLIB=../mdlib
HEADERS=$(MDLIB)/pvector.H $(MDLIB)/pmatrix.H $(MDLIB)/rpoly.H $(MDLIB)/boxes.H 
LDFLAGS=-lm -llapack -lblas #-g

all: mcsim veff

mcsim: mcsim_hardell.cpp $(HEADERS) $(MDLIB)/linked_cell_lists.H
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o mcsim_hardell mcsim_hardell.cpp 

veff: calc_veff_hardell.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o calc_veff_hardell calc_veff_hardell.cpp 
#gcc-9 -c ellipses.c -O3
#g++-9 -lm -llapack -lblas prova_ellipses.o  ellipses.o -o testhe
clean:
	rm -fr calc_veff_hardell mcsim_hardell
