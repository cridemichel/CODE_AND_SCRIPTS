#gauss Makefile
SDHDR=mono.h
LAMDIR=/usr/local/lam
UTIL = wtdutil heteroutil xvautil corutil jumps basins #basins_serial

gauss: link measures.o forces.o move.o init.o corutil basins cri2fra# $(UTIL)

link:
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/gauss/$(SDHDR) $(COMM)/mdsimdep.h 
# rigidABdbLJ utilities 
wtdutil: wtdutil.o
	$(CC) $(LDFLAGS) wtdutil.o -o $(BIN)/MOwtdutil -lm 

basins: basins.o forces.o $(COMM)/mdarray.o init.o
	mpicc $(LDFLAGS) basins.o $(COMM)/mdarray.o forces.o init.o -o $(BIN)/basins -L/usr/lib/lam/lib -lm 

cri2fra: cri2fra.o $(COMM)/mdarray.o
	$(CC) $(LDFLAGS) cri2fra.o $(COMM)/mdarray.o -o $(BIN)/cri2fraGauss -lm 

basins_serial: basins_serial.o forces.o $(COMM)/mdarray.o init.o
	mpicc $(LDFLAGS) basins_serial.o $(COMM)/mdarray.o forces.o init.o -o $(BIN)/basins -lm 

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

forces.o: forces.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o forces.o forces.c

cri2fra.o: $(COMM)/mdsimul.h $(COMM)/mdarray.c cri2fra.c
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o cri2fra.o cri2fra.c

basins.o: basins.c $(COMM)/mdarray.c forces.c $(SDHDR) $(COMM)/mdsimul.h init.c
	mpicc $(CFLAGS) $(MACROS) -c -I/usr/include/lam $(INCLUDES) -o basins.o basins.c

basins_serial.o: basins_serial.c $(COMM)/mdarray.c forces.c $(SDHDR) $(COMM)/mdsimul.h init.c
	gcc $(CFLAGS) $(MACROS) -c $(INCLUDES) -o basins_serial.o basins_serial.c
 
init.o: init.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o init.o init.c

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

