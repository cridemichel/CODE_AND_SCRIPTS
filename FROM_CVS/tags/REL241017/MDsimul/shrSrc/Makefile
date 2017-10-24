PERC = 	$(HOME)/shared
CC = gcc
CFLAGS = -O     

library:  ShrMem.h ShrMem.c
	$(CC) $(CFLAGS) -I$(PERC)/src -c -o $(PERC)/lib/libshare.a ShrMem.c

clean: 
	rm -f $(PERC)/lib/libshare.a








