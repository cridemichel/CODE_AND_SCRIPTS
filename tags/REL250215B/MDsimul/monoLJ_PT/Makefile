#monoLJ Makefile
SDHDR=mono_PT.h

UTIL = wtdutil heteroutil xvautil jumps corutil

monoLJ_PT: link measures.o forces_PT.o move_PT.o init_PT.o # $(UTIL)

link:
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/monoLJ_PT/$(SDHDR) $(COMM)/mdsimdep.h 
# rigidABdbLJ utilities 
wtdutil: wtdutil.o
	$(CC) $(LDFLAGS) wtdutil.o -o $(BIN)/MOwtdutil -lm 

heteroutil: heteroutil.o
	$(CC) $(LDFLAGS) heteroutil.o -o $(BIN)/MOheteroutil -lm 

xvautil: xvautil.o
	$(CC) $(LDFLAGS) xvautil.o -o $(BIN)/MOxvautil -lm 

corutil: corutil.o
	$(CC) $(LDFLAGS) corutil.o $(COMM)/scalProc_PT.o -o $(BIN)/MOcorutil -lm 

jumps: jump.o
	$(CC) $(LDFLAGS) jumps.o -o $(BIN)/MOjumps -lm 

measures.o: measures.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o measures.o measures.c

move_PT.o: move_PT.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o move_PT.o move_PT.c

forces_PT.o: forces_PT.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o forces_PT.o forces_PT.c

init_PT.o: init_PT.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o init_PT.o init_PT.c

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

