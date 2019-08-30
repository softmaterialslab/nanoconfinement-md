# This is a makefile.
# Use option -p in CC for profiling with gprof

PROG = md_simulation_confined_ions
BIN = ../bin/
bindir = bin
OBJ = main.o interface.o NanoconfinementMd.o functions.o md.o mdforces.o mdenergies.o
CCF = CC

# nanoHUB flags.
nanoHUBCC = mpicxx -O3 -g -Wall -fopenmp -std=c++11
nanoHUBLFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_program_options -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
nanoHUBCFLAG = -I${BOOST_INCDIR}
nanoHUBOFLAG = -o
# BigRed2 flags.
BigRed2CC = CC -O3 -g -Wall -fopenmp -std=c++11
BigRed2LFLAG = -lgsl -lgslcblas -lm -lboost_program_options -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
BigRed2CFLAG = -c
BigRed2OFLAG = -o
# General purpose flags.
CC = mpicxx -O3 -g -Wall -fopenmp -std=c++11
LFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_program_options -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
CFLAG = -c
OFLAG = -o

all: $(PROG)

install: create-dirs
	. /etc/environ.sh; use -e -r boost-1.62.0-mpich2-1.3-gnu-4.7.2; make CCF=nanoHUB all
	@echo "Installing $(PROG) into $(bindir) directory on nanohub"; cp -f $(PROG) $(BIN)

local-install: create-dirs
	make CCF=LOCAL all
	@echo "Installing $(PROG) into $(bindir) directory"; cp -f $(PROG) $(BIN)

cluster-install: create-dirs
	module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl; make CCF=BigRed2 all
	@echo "Installing $(PROG) into $(bindir) directory on the computing cluster"; cp -f $(PROG) $(BIN)

create-dirs:
	@echo "Checking and creating needed sub-directories in the $(BIN) directory"
	if ! [ -d $(BIN) ]; then mkdir $(BIN); fi
	if ! [ -d $(BIN)/outfiles ]; then mkdir $(BIN)/outfiles; fi
	if ! [ -d $(BIN)/infiles ]; then mkdir $(BIN)/infiles; fi
	if ! [ -d $(BIN)/data ]; then mkdir $(BIN)/data; fi
	@echo "Directory creation is over."

$(PROG) : $(OBJ)
ifeq ($(CCF),BigRed2)
	$(BigRed2CC) $(BigRed2OFLAG) $(PROG) $(OBJ) $(BigRed2LFLAG)
%.o : %.cpp
	$(BigRed2CC) -c $(BigRed2CFLAG) $< -o $@
else ifeq ($(CCF),nanoHUB)
	$(nanoHUBCC) $(nanoHUBOFLAG) $(PROG) $(OBJ) $(nanoHUBLFLAG)
%.o : %.cpp
	$(nanoHUBCC) -c $(nanoHUBCFLAG) $< -o $@
else
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LFLAG)
%.o : %.cpp
	$(CC) -c $(CFLAG) $< -o $@
endif

clean: dataclean
	rm -f *.o
	rm -f $(PROG)
	rm -rf $(BIN)/outfiles/*
	rm -rf $(BIN)/data/*
	rm -f $(BIN)/*.lammps
	rm -f $(BIN)/temp

dataclean:
	rm -f $(BIN)/outfiles/*.dat $(BIN)/outfiles/*.xyz  $(BIN)/outfiles/*.lammpstrj  $(BIN)/temp/*
	rm -f $(BIN)/data/*.dat $(BIN)/data/*.xyz  $(BIN)/data/*.lammpstrj
	rm -f $(BIN)/*.log
	rm -f $(BIN)/*.pbs
	rm -f $(BIN)/*.lammps

distclean: clean
	rm -f $(BIN)/$(PROG)

.PHONY: all install nanoHUB-install cluster-install distclean dataclean clean
