#bimix Makefile
SDHDR=bimixLMC.h

UTIL = wtdutil heteroutil bixva

bimix:  measuresLMC.o forcesLMC.o moveLMC.o initLMC.o md2asciiLMC link corutil cri2fra#$(UTIL)

link:
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/bimixLMC/$(SDHDR) $(COMM)/mdsimdep.h

corutil: corutil.o scalProc.o
	$(CC) $(LDFLAGS) -lm -o $(BIN)/bicorutilLMC scalProc.o corutil.o	

cri2fra: cri2fra.o $(COMM)/mdarrayMIX.o 
	$(CC) $(LDFLAGS) -lm cri2fra.o $(COMM)/mdarrayMIX.o -o $(BIN)/cri2fraMIX

md2asciiLMC: $(COMM)/md2ascii.o
	$(CC) $(LDFLAGS) $(INCLUDES) -lm -o $(BIN)/md2asciiLMC $(COMM)/md2ascii.o 

wtdutil: wtdutil.o
	$(CC) $(LDFLAGS) -lm wtdutil.o -o $(BIN)/wtdutil	

heteroutil: heteroutil.o
	$(CC) $(LDFLAGS) -lm heteroutil.o -o $(BIN)/heteroutil	

bixva: bixva.o
	$(CC) $(LDFLAGS) -lm bixva.o -o $(BIN)/bixva	

jump: jump.o
	$(CC) $(LDFLAGS) -lm jump.o -o $(BIN)/jump	


measuresLMC.o: measuresLMC.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o measuresLMC.o measuresLMC.c

moveLMC.o: moveLMC.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o moveLMC.o moveLMC.c

forcesLMC.o: forcesLMC.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o forcesLMC.o forcesLMC.c

initLMC.o: initLMC.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o initLMC.o initLMC.c

corutil.o: $(COMM)/mdsimul.h corutil.c
	 $(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o corutil.o corutil.c

cri2fra.o: $(COMM)/mdsimul.h cri2fra.c
	 $(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o cri2fra.o cri2fra.c

scalProc.o: scalProc.c $(COMM)/mdsimul.h $(COMM)/mdsimdep.h scalProc.h
	$(CC) $(CFLAGS) $(MACROS)  $(INCLUDES) -c  -o scalProc.o scalProc.c

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

