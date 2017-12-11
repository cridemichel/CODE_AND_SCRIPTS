#bimix Makefile
SDHDR=bimix.h
SDHDR1=bimix1.h
SDHDR2=bimix2.h
SDHDR3=bimix3.h
SDHDR4=bimix4.h

UTIL = wtdutil heteroutil bixva 

bimix:  md2asciiBM measures.o forces.o move.o init.o conjgrad.o mediafsk calcdrsq calcpress calcSk calcSkAB calcSkAA calcSkBB calcSkAAallk cri2fra fsk fk basinsSerial
# corutil basins$(UTIL)

link: 
	#simuln.sh $(PERC)/bimix/$(SDHDR) 0 
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/bimix/$(SDHDR) $(COMM)/mdsimdep.h

link1: 
	simuln.sh $(PERC)/bimix/$(SDHDR) 0 
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/bimix/$(SDHDR) $(COMM)/mdsimdep.h

link2: 
	simuln.sh $(PERC)/bimix/$(SDHDR) 2 
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/bimix/$(SDHDR) $(COMM)/mdsimdep.h

link3: 
	simuln.sh $(PERC)/bimix/$(SDHDR)  3
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/bimix/$(SDHDR) $(COMM)/mdsimdep.h
link4: 
	simuln.sh $(PERC)/bimix/$(SDHDR)  4
	rm $(COMM)/mdsimdep.h ; ln -s $(PERC)/bimix/$(SDHDR) $(COMM)/mdsimdep.h

basins: basins.o forces.o $(COMM)/mdarrayMIX.o init.o
	mpicc $(INCLUDES) $(LDFLAGS)  basins.o $(COMM)/mdarrayMIX.o forces.o init.o -o $(BIN)/basins -lm 
# mpicc_mpich $(LDFLAGS) -lm basins.o $(COMM)/mdarrayMIX.o forces.o init.o -o $(BIN)/basins
basinsSerial: basinsSerial.o forces.o $(COMM)/mdarrayMIX.o init.o
	$(CC) $(INCLUDES) $(LDFLAGS) basinsSerial.o $(COMM)/mdarrayMIX.o forces.o init.o -o $(BIN)/basinsSerial -lm

fsk: Fsk.o $(COMM)/mdarrayMIX.o init.o $(COMM)/mdio.o
	$(CC) $(INCLUDES) $(LDFLAGS) Fsk.o $(COMM)/mdarrayMIX.o init.o -o $(BIN)/buildfsk -lm

fk: Fk.o $(COMM)/mdarrayMIX.o init.o $(COMM)/mdio.o
	$(CC) $(INCLUDES) $(LDFLAGS) Fk.o $(COMM)/mdarrayMIX.o init.o -o $(BIN)/buildfk -lm

mediafsk: mediafsk.o 
	$(CC) $(LDFLAGS) -o $(BIN)/mediafsk mediafsk.o -lm

calcdrsq: calcdrsq.o forces.o $(COMM)/mdarrayMIX.o init.o 
	$(CC) $(LDFLAGS) calcdrsq.o init.o forces.o $(COMM)/mdarrayMIX.o -o $(BIN)/calcdrsq -lm 

calcpress: calcpress.o forces.o $(COMM)/mdarrayMIX.o init.o 
	$(CC) -g $(LDFLAGS)  calcpress.o init.o forces.o $(COMM)/mdarrayMIX.o -o $(BIN)/calcpress -lm 
calcSk: calcSk.o forces.o $(COMM)/mdarrayMIX.o init.o 
	$(CC) $(LDFLAGS)  calcSk.o init.o forces.o $(COMM)/mdarrayMIX.o -o $(BIN)/calcSk -lm 
calcSkAB: calcSkAB.o forces.o $(COMM)/mdarrayMIX.o init.o 
	$(CC) $(LDFLAGS)  calcSkAB.o init.o forces.o $(COMM)/mdarrayMIX.o -o $(BIN)/calcSkAB -lm 
calcSkBB: calcSkBB.o forces.o $(COMM)/mdarrayMIX.o init.o 
	$(CC) $(LDFLAGS) calcSkBB.o init.o forces.o $(COMM)/mdarrayMIX.o -o $(BIN)/calcSkBB  -lm 
calcSkAA: calcSkAA.o forces.o $(COMM)/mdarrayMIX.o init.o 
	$(CC) $(LDFLAGS) calcSkAA.o init.o forces.o $(COMM)/mdarrayMIX.o -o $(BIN)/calcSkAA -lm 
calcSkAAallk: calcSkAAallk.o forces.o $(COMM)/mdarrayMIX.o init.o 
	$(CC) $(LDFLAGS) calcSkAAallk.o init.o forces.o $(COMM)/mdarrayMIX.o -o $(BIN)/calcSkAAallk  -lm 


md2asciiBM: $(COMM)/md2ascii.o
	$(CC) $(LDFLAGS) $(INCLUDES)  -o $(BIN)/md2asciiBM $(COMM)/md2ascii.o -lm

corutil: corutil.o scalProc.o
	$(CC) $(LDFLAGS) -o $(BIN)/bicorutil scalProc.o corutil.o -lm 

cri2fra: cri2fra.o $(COMM)/mdarrayMIX.o 
	$(CC) $(LDFLAGS) cri2fra.o $(COMM)/mdarrayMIX.o -o $(BIN)/cri2fraMIX -lm 

wtdutil: wtdutil.o
	$(CC) $(LDFLAGS) wtdutil.o -o $(BIN)/wtdutil -lm 

heteroutil: heteroutil.o
	$(CC) $(LDFLAGS) heteroutil.o -o $(BIN)/heteroutil -lm 

bixva: bixva.o
	$(CC) $(LDFLAGS) bixva.o -o $(BIN)/bixva -lm 

jump: jump.o
	$(CC) $(LDFLAGS) jump.o -o $(BIN)/jump 	-lm 

measures.o: measures.c $(SDHDR) $(COMM)/mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o measures.o measures.c

basins.o: basins.c $(COMM)/mdarrayMIX.c forces.c $(SDHDR) $(COMM)/mdsimul.h init.c
	#mpicc_mpich $(CFLAGS) $(MACROS) -c $(INCLUDES) -o basins.o basins.c
	mpicc $(CFLAGS) $(MACROS) -c $(INCLUDES) -o basins.o basins.c

basinsSerial.o: basinsSerial.c $(COMM)/mdarrayMIX.c forces.c $(SDHDR) $(COMM)/mdsimul.h init.c
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o basinsSerial.o basinsSerial.c

calcdrsq.o: calcdrsq.c $(COMM)/mdarrayMIX.c init.c $(SDHDR) $(COMM)/mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o calcdrsq.o calcdrsq.c

calcpress.o: calcpress.c $(COMM)/mdarrayMIX.c init.c $(SDHDR) $(COMM)/mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS)   -c $(INCLUDES) -o calcpress.o calcpress.c
calcSk.o: calcSk.c $(COMM)/mdarrayMIX.c init.c $(SDHDR) $(COMM)/mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS)   -c $(INCLUDES) -o calcSk.o calcSk.c

calcSkAB.o: calcSkAB.c $(COMM)/mdarrayMIX.c init.c $(SDHDR) $(COMM)/mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS)   -c $(INCLUDES) -o calcSkAB.o calcSkAB.c
calcSkAA.o: calcSkAA.c $(COMM)/mdarrayMIX.c init.c $(SDHDR) $(COMM)/mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS)   -c $(INCLUDES) -o calcSkAA.o calcSkAA.c
calcSkBB.o: calcSkBB.c $(COMM)/mdarrayMIX.c init.c $(SDHDR) $(COMM)/mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS)   -c $(INCLUDES) -o calcSkBB.o calcSkBB.c
calcSkAAallk.o: calcSkAAallk.c $(COMM)/mdarrayMIX.c init.c $(SDHDR) $(COMM)/mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS)   -c $(INCLUDES) -o calcSkAAallk.o calcSkAAallk.c


Fsk.o: Fsk.c $(COMM)/mdarrayMIX.c $(SDHDR) $(COMM)/mdsimul.h init.c $(COMM)/mdio.c
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o Fsk.o Fsk.c

Fk.o: Fk.c $(COMM)/mdarrayMIX.c $(SDHDR) $(COMM)/mdsimul.h init.c $(COMM)/mdio.c
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o Fk.o Fk.c
	
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

#cri2fra.o: $(COMM)/mdsimul.h cri2fra.c
#	 $(CC) $(CFLAGS) $(MACROS) $(INCLUDES) -c -o cri2fra.o cri2fra.c

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

