#monoLJ Makefile
SDHDR=mono.h

UTIL = corutil

monoLJ: link measures.o forces.o move.o init.o corutil# $(UTIL)

link:
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/monoMC/$(SDHDR) $(COMM)/mdsimdep.h 
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

forces.o: forces.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o forces.o forces.c

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

