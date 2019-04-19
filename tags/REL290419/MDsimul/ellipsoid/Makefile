#ellipsoid Makefile
SDHDR=ellipsoid.h

UTIL = wtdutil heteroutil xvautil corutil jumps

ellipsoid: link measures.o move.o init.o corutil tree.o matrixalg.o nnl.o hesticky.o superellipsoid.o multiple_ll.o montecarlo.o hard_cylinders.o polysolv.o hard_ellipsoid.o
#parallel_ed.o # $(UTIL)

link:
	ln -sf $(PERC)/ellipsoid/$(SDHDR) $(COMM)/mdsimdep.h 
# rigidABdbLJ utilities 
wtdutil: wtdutil.o
	$(CC) $(LDFLAGS) wtdutil.o -o $(BIN)/MOwtdutil -lm 

heteroutil: heteroutil.o
	$(CC) $(LDFLAGS) heteroutil.o -o $(BIN)/MOheteroutil -lm 

xvautil: xvautil.o
	$(CC) $(LDFLAGS) xvautil.o -o $(BIN)/MOxvautil -lm 

corutil: corutil.o
	$(CC) $(LDFLAGS) corutil.o $(COMM)/scalProc.o -o $(BIN)/MOcorutil -lm 

jumps: jump.o
	$(CC) $(LDFLAGS) jumps.o -o $(BIN)/MOjumps -lm 

measures.o: measures.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o measures.o measures.c

move.o: move.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o move.o move.c

init.o: init.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o init.o init.c

tree.o: tree.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o tree.o tree.c
	
matrixalg.o: matrixalg.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o matrixalg.o matrixalg.c

nnl.o: nnl.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o nnl.o nnl.c
	
polysolv.o: polysolv.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o polysolv.o polysolv.c

hard_ellipsoid.o: hard_ellipsoid.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o hard_ellipsoid.o hard_ellipsoid.c
	
hesticky.o: hesticky.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o hesticky.o hesticky.c
	
superellipsoid.o: superellipsoid.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o superellipsoid.o superellipsoid.c

montecarlo.o: montecarlo.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o montecarlo.o montecarlo.c 
	
hard_cylinders.o: hard_cylinders.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o hard_cylinders.o hard_cylinders.c 
	
multiple_ll.o: multiple_ll.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o multiple_ll.o multiple_ll.c

parallel_dd.o: parallel_dd.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o parallel_dd.o parallel_dd.c

wtdutil.o: wtd.c wtd.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o wtdutil.o wtd.c

heteroutil.o: hetero.c hetero.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o hetero.o hetero.c

xvautil.o: xvautil.c xvautil.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o xvautil.o xvautil.c

jumps.o: jumps.c jumps.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o jumps.o jumps.c

corutil.o: corutil.c $(COMM)/mdsimul.h $(SDHDR) 
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o corutil.o corutil.c

realclean: clean

clean:
	rm -f *.o

