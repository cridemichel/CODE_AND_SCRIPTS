#bimix Makefile
SDHDR=monoLMC.h

UTIL =

monoLMC:  md2asciiMLMC link corutil measuresMLMC.o forcesMLMC.o moveMLMC.o initMLMC.o#$(UTIL)

link:
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/monoLMC/$(SDHDR) $(COMM)/mdsimdep.h

corutil: corutil.o scalProc.o
	$(CC) $(LDFLAGS) -lm -o $(BIN)/corutilMLMC scalProc.o corutil.o	

md2asciiMLMC: $(COMM)/md2ascii.o
	$(CC) $(LDFLAGS) $(INCLUDES) -lm -o $(BIN)/md2asciiLMC $(COMM)/md2ascii.o 


measuresMLMC.o: measuresMLMC.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o measuresMLMC.o measuresMLMC.c

moveMLMC.o: moveMLMC.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o moveMLMC.o moveMLMC.c

forcesMLMC.o: forcesMLMC.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o forcesMLMC.o forcesMLMC.c

initMLMC.o: initMLMC.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o initMLMC.o initMLMC.c

corutil.o: corutil.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o corutil.o corutil.c

scalProc.o: scalProc.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o scalProc.o scalProc.c

realclean: clean

clean:
	rm -f *.o

