#bimix Makefile
SDHDR=bimix.h
UTIL = wtdutil heteroutil bixva 

bimix:  md2asciiBM link measures.o forces.o move.o init.o conjgrad.o corutil cri2fra fsk mediafsk basins #$(UTIL)

link:
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/bimix/$(SDHDR) $(COMM)/mdsimdep.h

basins: basins.o forces.o $(COMM)/mdarrayMIX.o init.o
	mpicc $(INCLUDES )$(LDFLAGS) -lm basins.o $(COMM)/mdarrayMIX.o forces.o init.o -o $(BIN)/basins
# mpicc_mpich $(LDFLAGS) -lm basins.o $(COMM)/mdarrayMIX.o forces.o init.o -o $(BIN)/basins

fsk: Fsk.o $(COMM)/mdarrayMIX.o init.o $(COMM)/mdio.o
	$(CC) $(INCLUDES) $(LDFLAGS) -lm Fsk.o $(COMM)/mdarrayMIX.o init.o -o $(BIN)/buildfsk

mediafsk: mediafsk.o 
	$(CC) -lm $(LDFLAGS) -o mediafsk mediafsk.o

md2asciiBM: $(COMM)/md2ascii.o
	$(CC) $(LDFLAGS) $(INCLUDES) -lm -o $(BIN)/md2asciiBM $(COMM)/md2ascii.o 

corutil: corutil.o scalProc.o
	$(CC) $(LDFLAGS) -lm -o $(BIN)/bicorutil scalProc.o corutil.o	

cri2fra: cri2fra.o $(COMM)/mdarrayMIX.o 
	$(CC) $(LDFLAGS) -lm cri2fra.o $(COMM)/mdarrayMIX.o -o $(BIN)/cri2fraMIX

wtdutil: wtdutil.o
	$(CC) $(LDFLAGS) -lm wtdutil.o -o $(BIN)/wtdutil	

heteroutil: heteroutil.o
	$(CC) $(LDFLAGS) -lm heteroutil.o -o $(BIN)/heteroutil	

bixva: bixva.o
	$(CC) $(LDFLAGS) -lm bixva.o -o $(BIN)/bixva	

jump: jump.o
	$(CC) $(LDFLAGS) -lm jump.o -o $(BIN)/jump	

measures.o: measures.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o measures.o measures.c

basins.o: basins.c $(COMM)/mdarrayMIX.c forces.c $(SDHDR) $(COMM)/mdsimul.h init.c
	#mpicc_mpich $(CFLAGS) $(MACROS) -c $(INCLUDES) -o basins.o basins.c
	mpicc $(CFLAGS) $(MACROS) -c $(INCLUDES) -o basins.o basins.c

Fsk.o: Fsk.c $(COMM)/mdarrayMIX.c $(SDHDR) $(COMM)/mdsimul.h init.c $(COMM)/mdio.c
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o Fsk.o Fsk.c
	
mediafsk.o: mediafsk.c 
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mediafsk.o mediafsk.c

move.o: move.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o move.o move.c

forces.o: forces.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o forces.o forces.c

init.o: init.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o init.o init.c

conjgrad.o: conjgrad.c $(COMM)/mdsimul.h $(SDHDR)
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o conjgrad.o conjgrad.c

corutil.o: $(COMM)/mdsimul.h corutil.c
	 $(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o corutil.o corutil.c

cri2fra.o: $(COMM)/mdsimul.h cri2fra.c
	 $(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o cri2fra.o cri2fra.c

scalProc.o: scalProc.c $(COMM)/mdsimul.h $(COMM)/mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS)  $(INCLUDES) -c  -o scalProc.o scalProc.c

cri2fra.o: cri2fra.c $(COMM)/mdarrayMIX.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o cri2fra.o cri2fra.c

wtdutil.o: wtdconstr.c wtdconstr.h $(SDHDR) $(COMM)/mdsimul.h
	echo -n "eccomi qua" $(INCLUDES)
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o wtdutil.o wtdconstr.c

heteroutil.o: heteroconstr.c heteroconstr.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o heteroutil.o heteroconstr.c

bixva.o: bixva.c bixva.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o bixva.o bixva.c

jump.o: jumpconstr.c jumpconstr.h $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o jump.o jumpconstr.c


realclean: clean

clean:
	rm -f *.o

