#rigidABdbLJ Makefile
SDHDR=laponite.h

UTIL = #xvautil wtdutil heteroutil

laponite:  link measures.o forces.o move.o init.o $(UTIL)

link:
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/laponite/$(SDHDR) $(COMM)/mdsimdep.h 
# rigidABdbLJ utilities 
wtdutil: wtdutil.o
	$(CC) $(LDFLAGS) wtdutil.o -o $(BIN)/wtdutil -lm 

heteroutil: heteroutil.o
	$(CC) $(LDFLAGS) heteroutil.o -o $(BIN)/heteroutil -lm 

xvautil: xvautil.o
	$(CC) $(LDFLAGS) xvautil.o -o $(BIN)/xvautil -lm 

jump: jump.o
	$(CC) $(LDFLAGS) jump.o -o $(BIN)/jump	-lm 

measures.o: measures.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o measures.o measures.c

move.o: move.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o move.o move.c

forces.o: forces.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o forces.o forces.c

init.o: init.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o init.o init.c

wtdutil.o: wtdconstr.c wtdconstr.h $(SDHDR) $(COMM)/mdsimul.h
#	echo -n "eccomi qua" $(INCLUDES)
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o wtdutil.o wtdconstr.c

heteroutil.o: heteroconstr.c heteroconstr.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o heteroutil.o heteroconstr.c

xvautil.o: xvautil.c xvautil.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o xvautil.o xvautil.c

jump.o: jumpconstr.c jumpconstr.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o jump.o jumpconstr.c


realclean: clean

clean:
	rm -f *.o

