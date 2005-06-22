CC = gcc
#CFLAGS = -Wall -g -pipe -D_DEBUG
#CFLAGS = -Wall -pipe
#CFLAGS = -Wall -pedantic
#CFLAGS = -Wall 
CFLAGS = 
ARCH =
LIB  = -lm
OPT  = -g
#OPT  = -O2

# headers and routine to read cristiano's ellipsoids files
Hread_crist = struct_conf_cristiano.h
Pread_crist = io_cristiano_ellipses.c

#headers and routines for the spherical harmonics
Hspher_harm= 
Pspher_harm= 

#headers and routines for rotation matrices
Hrotmatrix=
Protmatrix=


# defines for the R2 (MSD) program
DEF_R2 = -DNread=195  -DNpart=256
#DEF_R2 = -DNread=195 -DNpart=256 -DSCHILLING=1 

# defines for the SQY program 
DEF_SQY = -DNpart=256

# defines for the CORR-Y program 
DEF_CORR_Y = -DNpart=256

# defines for the CORR-SQY program 
DEF_CORR_SQY = -DNpart=256

#defines for the ORDPAR program
DEF_SHOW = -DNpart=256

#defines for the SHOW program
DEF_SHOW = -DNpart=256

#defines for the CLUSTER program
DEF_CLUSTER = -DNpart=256

#defines for the CLUSTER program
DEF_COUPLES = -DNpart=256

default:
	@echo
	@echo " ****** This Makefile has the following options:"
	@echo
	@echo "   make R2      : calculates the mean square displacement"
	@echo "   make SQY     :(fourier-transforms * Ylm) of a configuration" 
	@echo "   make AVE-SQY : averages the files calculated with SQY.x"
	@echo "   make CORR-Y  : correlation functions for Y10 and Y20"
	@echo "   make CORR-SQY: correlation functions for SQY at a given Q-vector"
	@echo "   make ORDPAR  : nematic order parameter"
	@echo "   make CLUSTER : processe a configuration for molgl" 
	@echo "   make SHOW    : show infos about a configuration" 
	@echo "   make clean   : cleans files"
	@echo 
	@echo "   for more details, read REM.TXT"
	@echo 
	@echo 

# 
#all: R2 SQY CORR-Y CORR-SQY CLUSTER
all: COUPLES


# calculate R2
R2: $(Pread_crist)  $(Hread_crist) calc_R2_cristiano.c
	$(CC) $(CFLAGS) $(ARCH) $(LIB) $(INCL) $(OPT) $(Pread_crist)  $(DEF_R2) calc_R2_cristiano.c -o R2.x


#calculate SQY
SQY: $(Pread_crist)  $(Hread_crist) calc_SQY_cristiano.c
	$(CC) $(CFLAGS) $(ARCH) $(LIB) $(INCL) $(OPT) $(Pread_crist) $(DEF_SQY) calc_SQY_cristiano.c -o SQY.x

#averages the files calculated with SQY.x
AVE-SQY:
	cc -lm ave_SQY.c -o AVE-SQY.x

# calculate correlation in Y10, Y20
CORR-Y: $(Pread_crist)  $(Hread_crist) calc_corr_Y_cristiano.c
	$(CC) $(CFLAGS) $(ARCH) $(LIB) $(INCL) $(OPT) $(Pread_crist)  $(DEF_CORR_Y) calc_corr_Y_cristiano.c -o CORR-Y.x

# calculate correlation in SQY for a given Q-vector
CORR-SQY: $(Pread_crist)  $(Hread_crist) calc_corr_SQY_cristiano.c
	$(CC) $(CFLAGS) $(ARCH) $(LIB) $(INCL) $(OPT) $(Pread_crist)  $(DEF_CORR_SQY) calc_corr_SQY_cristiano.c -o CORR-SQY.x


#nematic order parameter
ORDPAR: $(Pread_crist)  $(Hread_crist) $(Pspher_harm) $(Protmatrix) calc_ordpar_cristiano.c
	$(CC) $(CFLAGS) $(ARCH) $(LIB) $(INCL) $(OPT) $(Pread_crist)  $(Pspher_harm) $(Protmatrix) $(DEF_SHOW) calc_ordpar_cristiano.c -o ORDPAR.x

#show couples
COUPLES: $(Pread_crist)  $(Hread_crist) $(Pspher_harm) $(Protmatrix) show_couples_cristiano.c
	$(CC) $(CFLAGS) $(ARCH) $(LIB) $(INCL) $(OPT) $(Pread_crist)  $(Pspher_harm) $(Protmatrix) $(DEF_CLUSTER) show_couples_cristiano.c -o COUPLES.x

#show clusters
CLUSTER: $(Pread_crist)  $(Hread_crist) $(Pspher_harm) $(Protmatrix) show_cluster_cristiano.c
	$(CC) $(CFLAGS) $(ARCH) $(LIB) $(INCL) $(OPT) $(Pread_crist)  $(Pspher_harm) $(Protmatrix) $(DEF_CLUSTER) show_cluster_cristiano.c -o CLUSTER.x

#show infos
SHOW: $(Pread_crist)  $(Hread_crist) $(Pspher_harm) $(Protmatrix) show_params_cristiano.c
	$(CC) $(CFLAGS) $(ARCH) $(LIB) $(INCL) $(OPT) $(Pread_crist)  $(Pspher_harm) $(Protmatrix) $(DEF_SHOW) show_params_cristiano.c -o SHOW.x


clean:
	rm -fv *.o a.out
	rm -fv *~* \#*
	rm -fv *.x
