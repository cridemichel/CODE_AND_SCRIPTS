#MAIN Makefile
# === >>> Put your PERC here !!! <<< ===
export PERC=$(HOME)/MDsimul
# ======================================
export COMM = $(PERC)/commSrc

export CC=gcc
#export CC=mpicc
export MPIDIR =/usr/include/lam
#export CC=mpcc_r#per SP4
#export MPIDIR=$(MPICH)#per beowulf 
#export CC= $(MPIDIR)/bin/hcc
#export CC=mpicc_mpich
#export CFLAGS  = -Wall -mcpu=i686 -O3\
	  -frerun-cse-after-loop -fstrength-reduce -fexpensive-optimizations\
          -felide-constructors -malign-loops=2 -malign-jumps=2\
	  -malign-functions
export CFLAGS  = -Wall  -O3  -mcpu=athlon-xp -march=athlon-xp 
#export CFLAGS  =-Wall -O3 -msse2 -mfpmath=sse -march=pentium4 -mcpu=pentium4 
#export CFLAGS  =-Wall -O3  -march=pentium4 -mcpu=pentium4 
#export CFLAGS=  -B -tp p6 -fast#ottimizzazioni per beowulf
#export CFLAGS= -Wall -O3 -mcpu=i686 -march=i686
#export CFLAGS=  -Wall -O3 -mcpu=athlon-xp -march=athlon-xp
#CFLAGS = -Wall -m486 -fomit-frame-pointer -finline-functions -O2
#export CFLAGS  = -qarch=pwr4 -qtune=pwr4  -O3 -bmaxdata:0x70000000 -qcpluscmt #-I$(TURBOMP)/POWER4/include#per SP4
#export CFLAGS= -g -Wall 
# MACROS passed to the pre-processor
export MACROS= -UMD_SQWELL -UMD_HSVISCO -UMD_RESPA_FILTER -UMD_RAPACONSTR -UMD_RESPA_NOSELONG -UMD_RESPA_NPT -UMD_RESPA_SWITCH -UMD_FENE -UMD_RESPA -DMD_LOADMESH -DXDISP -UATPTENS -DMOLPTENS -UATPRESS #-DMD_INFBARRIER#-DMD_SQWELL# -DSOFT_SPHERE  -DMD_STATIC_PP8 # -DMD_BONDHIST -DMD_BILOG -DMPI##-DMD_TSHAPED# -DSOFT_SPHERE  -DMD_STATIC_PP6#-DMDLLINT#-DMPI##-DMC_SIMUL -DSOFT_SPHERE #-DMPI -DMD_GRAVITY -DPRINTER_PROC -DTAPPING -DMD_LOADMESH
#export MACROS= -DMD_LOADMESH -DSOFT_SPHERE -DMD_STATIC_PP18 -DXDISP  -DATPTENS -UMOLPTENS -UATPRESS -DORIGIN -DMPI #-DPRINTER_PROC -DTAPPING -DMD_LOADMESH #per SP4
export LDFLAGS=#-L /usr/lib/lam/lib -llamf77mpi -lmpi -llam -lnsl -lutil -lpmpi -llammpio
#export LDFLAGS= #-L$(TURBOMP)/POWER4/lib -lxlf90#per SP4
# ATPTENS defined(-DATPTENS -U... undefined)=>calculate atomic pressure tensor
# MOLPTENS defined => calculate mlecular pressor tensor
# ATPRESS defined => calculate atomic pressure(if ATPTENS defined the atomic
# pressure is automatically calculate)

export BIN = $(PERC)/bin

#directories to search for #include
export INCLUDES =  -I$(MPIDIR)  -I. -I$(PERC)/include -I$(COMM)

# ============================================================================
COMMOBJ = $(COMM)/mdsimul.o $(COMM)/mdinit.o $(COMM)/mdarray.o $(COMM)/mdutil.o $(COMM)/mdio.o $(COMM)/mdchoose.o

COMMOBJMIX = $(COMM)/mdsimul.o $(COMM)/mdinit.o $(COMM)/mdarrayMIX.o $(COMM)/mdutil.o $(COMM)/mdio.o $(COMM)/mdchoose.o 

COMMONMLMC = $(COMM)/mdsimul.o $(COMM)/mdinitMLMC.o $(COMM)/mdarrayMLMC.o $(COMM)/mdutil.o $(COMM)/mdio.o $(COMM)/mdchoose.o 

COMMOBJ_PT = $(COMM)/mdsimul_PT.o $(COMM)/mdinit_PT.o $(COMM)/mdarray_PT.o $(COMM)/mdutil.o $(COMM)/mdio.o $(COMM)/mdchoose.o

COMMOBJ_DPT = $(COMM)/mdsimul_PT.o $(COMM)/mdinit_DPT.o $(COMM)/mdarray_DPT.o $(COMM)/mdutil.o $(COMM)/mdio_DPT.o $(COMM)/mdchoose_DPT.o

LAPO=$(PERC)/laponite
LAPOOBJS= $(LAPO)/move.o $(LAPO)/measures.o $(LAPO)/init.o $(LAPO)/forces.o

NMPOLY=$(PERC)/nmpoly
NMPOLYOBJS= $(NMPOLY)/move.o $(NMPOLY)/measures.o $(NMPOLY)/init.o $(NMPOLY)/forces.o $(NMPOLY)/moveRespa.o


# dependencies for all simulation codes
# rigid AB dumbell interacting by LJ12-6 potential
RIGID = $(PERC)/rigidABdbLJ

RIGIDOBJS = $(RIGID)/move.o $(RIGID)/measures.o $(RIGID)/init.o $(RIGID)/forces.o

BIMIX = $(PERC)/bimix

BIMIXOBJS = $(BIMIX)/move.o $(BIMIX)/measures.o $(BIMIX)/init.o $(BIMIX)/forces.o $(BIMIX)/conjgrad.o

BIMIXLMC= $(PERC)/bimixLMC

BIMIXLMCOBJS = $(BIMIXLMC)/moveLMC.o $(BIMIXLMC)/measuresLMC.o $(BIMIXLMC)/initLMC.o $(BIMIXLMC)/forcesLMC.o

MONOLMC=$(PERC)/monoLMC
MONOLMCOBJS = $(MONOLMC)/moveMLMC.o $(MONOLMC)/measuresMLMC.o $(MONOLMC)/initMLMC.o $(MONOLMC)/forcesMLMC.o

#monoatomic Lennard-Jones system with Ruocco many body term
MONO = $(PERC)/monoLJ

MONOOBJS =  $(MONO)/move.o $(MONO)/measures.o $(MONO)/init.o $(MONO)/forces.o

#monoatomic Monte Carlo  
MONOMC = $(PERC)/monoMC

MONOMCOBJS =  $(MONOMC)/move.o $(MONOMC)/measures.o $(MONOMC)/init.o $(MONOMC)/forces.o


#gauss
GAUSS = $(PERC)/gauss

GAUSSOBJS =  $(GAUSS)/move.o $(GAUSS)/measures.o $(GAUSS)/init.o $(GAUSS)/forces.o
#monoatomic hard sphere
MONOHS= $(PERC)/monohs
MONOHSOBJS= $(MONOHS)/move.o $(MONOHS)/measures.o $(MONOHS)/init.o $(MONOHS)/tree.o
#biatomic mixture of hard spheres
BIMIXHS= $(PERC)/bimixhs
BIMIXHSOBJS= $(BIMIXHS)/move.o $(BIMIXHS)/measures.o $(BIMIXHS)/init.o $(BIMIXHS)/tree.o

#biatomic mixture of hard spheres
ELLIPSOID= $(PERC)/ellipsoid
ELLIPSOIDOBJS= $(ELLIPSOID)/move.o $(ELLIPSOID)/measures.o $(ELLIPSOID)/init.o $(ELLIPSOID)/tree.o

#monoatomic Lennard-Jones with Ruocco many body term and Parallel Tempering
MONO_PT = $(PERC)/monoLJ_PT

MONO_PT_OBJS =  $(MONO_PT)/move_PT.o $(MONO_PT)/measures.o $(MONO_PT)/init_PT.o $(MONO_PT)/forces_PT.o

#monoatomic Lennard-Jones with Ruocco many body term and Parallel Tempering
MONO_DPT = $(PERC)/monoLJ_DPT

MONO_DPT_OBJS =  $(MONO_DPT)/move_DPT.o $(MONO_DPT)/measures.o $(MONO_DPT)/init_DPT.o $(MONO_DPT)/forces_DPT.o

default: monoLJ

all: rigidABdbLJ monoLJ H2S monoLJPT monoLJDPT bimix bimixLMC monoLMC monoMC

# make sub-make
rigidABdbLJ: rigidAB_link common rigidABdbLJObjs#common always here
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ) $(RIGIDOBJS) -o $(BIN)/rigidABdbLJ -lm 

bimix: bimix_link commonMIX bimixObjs#common always here
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJMIX) $(BIMIXOBJS) -o $(BIN)/bimix -lm
	
bimix1:  bimix_link1 commonMIX bimixObjs#common always here
	$(CC) $(COMMOBJMIX) $(BIMIXOBJS) -o $(BIN)/bimix2 $(LDFLAGS) -L$(PERC)/lib/ -lm 
bimix2: bimix_link2 commonMIX bimixObjs#common always here
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJMIX) $(BIMIXOBJS) -o $(BIN)/bimix2 -lm 
bimix3: bimix_link3 commonMIX bimixObjs#common always here
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJMIX) $(BIMIXOBJS) -o $(BIN)/bimix3 -lm 
bimix4: bimix_link4 commonMIX bimixObjs#common always here
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJMIX) $(BIMIXOBJS) -o $(BIN)/bimix4 -lm 

bimixLMC: bimixLMC_link commonMIX bimixLMCObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJMIX) $(BIMIXLMCOBJS) -o $(BIN)/bimixLMC -lm 

laponite: lapo_link common lapoObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ) $(LAPOOBJS) -o $(BIN)/laponite -lm 

nmpoly: nmpoly_link common nmpolyObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ) $(NMPOLYOBJS) -o $(BIN)/nmpoly -lm 


monoLJ: monoLJ_link common monoLJObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ) $(MONOOBJS) -o $(BIN)/monoLJ -lm 

monoMC: monoMC_link common monoMCObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ) $(MONOMCOBJS) -o $(BIN)/monoMC -lm 


gauss: gauss_link common gaussObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ) $(GAUSSOBJS) -o $(BIN)/gauss -lm 

monohs: monohs_link common monohsObjs
	$(CC) -L$(PERC)/lib/ $(COMMOBJ) $(MONOHSOBJS) -o $(BIN)/monohs $(LDFLAGS) -lm 

bimixhs: bimixhs_link common bimixhsObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ) $(BIMIXHSOBJS) -o $(BIN)/bimixhs -lm 

ellipsoid: ellipsoid_link common ellipsoidObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ) $(ELLIPSOIDOBJS) -o $(BIN)/ellipsoid -lm 

monoLMC: monoLMC_link commonMLMC monoLMCObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMONMLMC) $(MONOLMCOBJS) -o $(BIN)/monoLMC -lm 

monoLJPT: monoLJPT_link commonPT monoLJPTObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ_PT) $(MONO_PT_OBJS) -o $(BIN)/monoLJPT -lm 

monoLJDPT: monoLJDPT_link commonDPT monoLJDPTObjs
	$(CC) $(LDFLAGS) -L$(PERC)/lib/ $(COMMOBJ_DPT) $(MONO_DPT_OBJS) -o $(BIN)/monoLJDPT -lm  

H2S: common

rigidABdbLJObjs:
	cd $(PERC)/rigidABdbLJ && $(MAKE) 

bimixObjs:
	cd $(PERC)/bimix && $(MAKE) 

bimixLMCObjs:
	cd $(PERC)/bimixLMC && $(MAKE) 

lapoObjs:
	cd $(PERC)/laponite && $(MAKE) 

nmpolyObjs:
	cd $(PERC)/nmpoly && $(MAKE) 

monoLMCObjs:
	cd $(PERC)/monoLMC && $(MAKE) 

gaussObjs:
	cd $(PERC)/gauss && $(MAKE) 


H2SObjs:
	cd $(PERC)/H2S && $(MAKE)

monoLJObjs:
	cd $(PERC)/monoLJ && $(MAKE)

monoMCObjs:
	cd $(PERC)/monoMC && $(MAKE)

monohsObjs:
	cd $(PERC)/monohs && $(MAKE)

bimixhsObjs:
	cd $(PERC)/bimixhs && $(MAKE)

ellipsoidObjs:
	cd $(PERC)/ellipsoid && $(MAKE)

monoLJPTObjs:
	cd $(PERC)/monoLJ_PT && $(MAKE)

monoLJDPTObjs:
	cd $(PERC)/monoLJ_DPT && $(MAKE)

common: 
	cd $(COMM) && $(MAKE)

commonMIX: 
	cd $(COMM) && $(MAKE) commonMIX

commonMLMC:
	cd $(COMM) && $(MAKE) commonMLMC

commonPT:
	cd $(COMM) && $(MAKE) commonPT

commonDPT: 
	cd $(COMM) && $(MAKE) commonDPT

# if you change the simulation code change here!!!!!!!!!!!!!!
link:
	cd $(PERC)/bimix && $(MAKE) link

bimix_link:
	cd $(PERC)/bimix && $(MAKE) link

bimix_link1:
	cd $(PERC)/bimix && $(MAKE) link1
bimix_link2:
	cd $(PERC)/bimix && $(MAKE) link2
bimix_link3:
	cd $(PERC)/bimix && $(MAKE) link3
bimix_link4:
	cd $(PERC)/bimix && $(MAKE) link4

bimixLMC_link:
	cd $(PERC)/bimixLMC && $(MAKE) link

lapo_link:
	cd $(PERC)/laponite && $(MAKE) link

nmpoly_link:
	cd $(PERC)/nmpoly && $(MAKE) link


monoLMC_link:
	cd $(PERC)/monoLMC && $(MAKE) link

monoLJDPT_link:
	cd $(PERC)/monoLJ_DPT && $(MAKE) link

monoLJPT_link:
	cd $(PERC)/monoLJ_PT && $(MAKE) link

monoLJ_link:
	cd $(PERC)/monoLJ && $(MAKE) link

monoMC_link:
	cd $(PERC)/monoMC && $(MAKE) link


monohs_link:
	cd $(PERC)/monohs && $(MAKE) link

bimixhs_link:
	cd $(PERC)/bimixhs && $(MAKE) link

ellipsoid_link:
	cd $(PERC)/ellipsoid && $(MAKE) link

gauss_link:
	cd $(PERC)/gauss && $(MAKE) link

rigidAB_link:
	cd $(PERC)/rigidABdbLJ  && $(MAKE) link

realclean: clean 

# delete all bins
realclean: clean
	rm -f $(BIN)/*

clean:
	cd rigidABdbLJ && $(MAKE) clean
	cd gauss && $(MAKE) clean
	cd bimix && $(MAKE) clean
	cd bimixLMC && $(MAKE) clean
	cd monoLJ && $(MAKE) clean
	cd monohs && $(MAKE) clean
	cd bimixhs && $(MAKE) clean
	cd commSrc && $(MAKE) clean
	cd monoLJ_PT && $(MAKE) clean
	cd monoLJ_DPT && $(MAKE) clean	
	cd monoLMC && $(MAKE) clean
	cd laponite && $(MAKE) clean
	cd nmpoly && $(MAKE) clean

