#commSrc Makefile
COMMEXECS = md2ascii bak2cor calcEta# log calcEta chkStatus 
COMMEXECS_PT = md2ascii_PT log calcEta chkStatus bak2cor
COMMEXECS_DPT = md2ascii_DPT calcEta chkStatus bak2cor xva2ascii
OBJS= mdutil.o mdsimul.o mdio.o mdchoose.o mdarray.o mdinit.o scalProc.o bak2cor.o
OBJSMIX= mdutil.o mdsimul.o mdio.o mdchoose.o mdarrayMIX.o mdinit.o bak2cor.o
OBJS_MLMC = mdutil.o mdsimul.o mdio.o mdchoose.o mdarrayMLMC.o mdinitMLMC.o bak2cor.o
OBJSPT= mdutil.o mdsimul_PT.o mdio.o mdchoose.o mdarray_PT.o mdinit_PT.o bak2cor.o#scalProc_PT.o
OBJSDPT= mdutil.o mdsimul_PT.o mdio_DPT.o mdchoose_DPT.o mdarray_DPT.o mdinit_DPT.o bak2cor.o#scalProc_PT.o

common: $(COMMEXECS) $(OBJS)
commonMIX: $(COMMEXECS) $(OBJSMIX)
commonPT: $(COMMEXECS_PT) $(OBJSPT)
commonDPT: $(COMMEXECS_DPT) $(OBJSDPT)
commonMLMC: $(COMMEXECS) $(OBJS_MLMC)

# generate executable directly
chkStatus: chkStatus.o mdio.o mdutil.o
	$(CC)  $(LDFLAGS)  $(INCLUDES) -o chkStatus chkStatus.o mdio.o mdutil.o

log: log.c $(COMM)/mdsimdep.h $(COMM)/mdsimul.h
	$(CC) $(LDFLAGS) $(CFLAGS) $(INCLUDES) log.c -o $(BIN)/log

md2ascii: md2ascii.o
	$(CC) $(LDFLAGS) $(INCLUDES)  -o $(BIN)/md2ascii md2ascii.o -lm

md2ascii_PT: md2ascii_PT.o
	$(CC) $(LDFLAGS) $(INCLUDES) -o $(BIN)/md2ascii md2ascii_PT.o -lm 

md2ascii_DPT: md2ascii_DPT.o
	$(CC) $(LDFLAGS) $(INCLUDES) -o $(BIN)/md2ascii md2ascii_DPT.o -lm 

calcEta: calcEta.o
	$(CC) $(LDFLAGS) calcEta.o -o $(BIN)/calcEta -lm	

xva2ascii: xva2ascii.o
	$(CC) $(LDFLAGS) xva2ascii.o -o $(BIN)/xva2ascii -lm		

bak2cor: bak2cor.o
	$(CC) $(LDFLAGS) bak2cor.o -o $(BIN)/bak2cor -lm	

mdio.o : mdio.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o mdio.o mdio.c 

mdio_DPT.o : mdio_DPT.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o mdio_DPT.o mdio_DPT.c 

mdchoose.o : mdchoose.c mdsimul.h mdsimdep.h 
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o mdchoose.o mdchoose.c 

mdchoose_DPT.o : mdchoose_DPT.c mdsimul.h mdsimdep.h 
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o mdchoose_DPT.o mdchoose_DPT.c 

mdsimul.o: mdsimul.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdsimul.o mdsimul.c

mdsimul_PT.o: mdsimul_PT.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdsimul_PT.o mdsimul_PT.c

md2ascii.o: md2ascii.c mdsimdep.h mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o md2ascii.o md2ascii.c 

md2ascii_PT.o: md2ascii_PT.c mdsimdep.h mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o md2ascii_PT.o md2ascii_PT.c 

md2ascii_DPT.o: md2ascii_DPT.c mdsimdep.h mdsimul.h 
	$(CC) $(CFLAGS) $(MACROS) $(INCLUDES)  -c -o md2ascii_DPT.o md2ascii_DPT.c 

scalProc.o: scalProc.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS)  $(INCLUDES) -c  -o scalProc.o scalProc.c

scalProc_PT.o: scalProc_PT.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS)  $(INCLUDES) -c  -o scalProc_PT.o scalProc_PT.c

mdinit.o: mdinit.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdinit.o mdinit.c

mdinitMLMC.o: mdinitMLMC.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdinitMLMC.o mdinitMLMC.c

mdinit_PT.o: mdinit_PT.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdinit_PT.o mdinit_PT.c

mdinit_DPT.o: mdinit_DPT.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdinit_DPT.o mdinit_DPT.c

mdarray.o : mdarray.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdarray.o mdarray.c

mdarrayMIX.o : mdarrayMIX.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdarrayMIX.o mdarrayMIX.c

mdarrayMLMC.o : mdarrayMLMC.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdarrayMLMC.o mdarrayMLMC.c

mdarray_PT.o : mdarray_PT.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdarray_PT.o mdarray_PT.c

mdarray_DPT.o : mdarray_DPT.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdarray_DPT.o mdarray_DPT.c

mdutil.o: mdutil.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o mdutil.o mdutil.c

bak2cor.o: bak2cor.c mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o bak2cor.o bak2cor.c

calcEta.o: calcEta.c calcEta.h mdsimul.h mdsimdep.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o calcEta.o calcEta.c

xva2ascii.o: xva2ascii.c mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o xva2ascii.o xva2ascii.c

chkStatus.o: chkStatus.c mdsimul.h
	$(CC) $(CFLAGS) $(MACROS) -c $(INCLUDES) -o chkStatus.o chkStatus.c

clean:
	rm -f *.o

realclean: 
	rm -f *.o 
