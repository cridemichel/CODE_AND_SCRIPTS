#monoLJ Makefile
SDHDR=mono_DPT.h

UTIL = wtdutil diffusion heteroutil xvautil jumps corutil basins buildPE

monoLJ_DPT: link measures.o forces_DPT.o move_DPT.o init_DPT.o diffusion basins #cri2fra $(UTIL)

link:
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/monoLJ_DPT/$(SDHDR) $(COMM)/mdsimdep.h 

# rigidABdbLJ utilities 
basins: basins.o forces_DPT.o $(COMM)/mdarray_DPT.o init_DPT.o
	$(CC) $(LDFLAGS) basins.o $(COMM)/mdarray_DPT.o forces_DPT.o init_DPT.o -o $(BIN)/basins -lm 

buildPE: buildPE.o forces_DPT.o $(COMM)/mdarray_DPT.o init_DPT.o
	gcc $(INCLUDES) $(LDFLAGS) buildPE.o $(COMM)/mdarray_DPT.o forces_DPT.o init_DPT.o -o $(BIN)/buildPE -lm 

basins_serial: basins_serial.o forces_DPT.o $(COMM)/mdarray_DPT.o init_DPT.o
	gcc $(LDFLAGS) basins_serial.o $(COMM)/mdarray_DPT.o forces_DPT.o init_DPT.o -o $(BIN)/basins -lm 

cri2fra: cri2fra.o $(COMM)/mdarray_DPT.o 
	$(CC) $(LDFLAGS) cri2fra.o $(COMM)/mdarray_DPT.o -o $(BIN)/cri2fra -lm 

diffusion: diffusion.o $(COMM)/mdarray_DPT.o init_DPT.o forces_DPT.o
	 gcc $(LDFLAGS) diffusion.o  init_DPT.o forces_DPT.o $(COMM)/mdarray_DPT.o -o $(BIN)/diffusion -lm 

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

move_DPT.o: move_DPT.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o move_DPT.o move_DPT.c

forces_DPT.o: forces_DPT.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o forces_DPT.o forces_DPT.c

init_DPT.o: init_DPT.c  $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o init_DPT.o init_DPT.c

basins.o: basins.c $(COMM)/mdarray_DPT.c forces_DPT.c $(SDHDR) $(COMM)/mdsimul.h init_DPT.c
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o basins.o basins.c

buildPE.o: buildPE.c $(COMM)/mdarray_DPT.c forces_DPT.c $(SDHDR) $(COMM)/mdsimul.h init_DPT.c
	gcc $(CFLAGS) $(MACROS) -c $(INCLUDES) -o buildPE.o buildPE.c

basins_serial.o: basins_serial.c $(COMM)/mdarray_DPT.c forces_DPT.c $(SDHDR) $(COMM)/mdsimul.h init_DPT.c
	gcc $(CFLAGS) $(MACROS) -c $(INCLUDES) -o basins_serial.o basins_serial.c
 
cri2fra.o: cri2fra.c $(COMM)/mdarray_DPT.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o cri2fra.o cri2fra.c

diffusion.o: diffusion.c diffusion.h $(COMM)/mdsimul.h
	gcc $(CFLAGS) $(MACROS) -c $(INCLUDES) -o diffusion.o diffusion.c

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

