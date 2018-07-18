# This is a makefile.
# Use option -p in CC for profiling with gprof

PROG = md_simulation_confined_ions
BIN = ../bin/
APPHOME = ../
bindir = bin
homedir = home
OBJ = main.o interface.o NanoconfinementMd.o functions.o md.o mdforces.o mdenergies.o
CCF = CC

# nanoHUB flags. 
nanoHUBCC = mpicxx -O3 -g -Wall -fopenmp
nanoHUBLFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_program_options -lboost_mpi -lboost_serialization
nanoHUBCFLAG = -I${BOOST_INCDIR}
nanoHUBOFLAG = -o
# BigRed2 flags. 
BigRed2CC = CC -O3 -g -Wall -fopenmp
BigRed2LFLAG = -lgsl -lgslcblas -lm -lboost_program_options -lboost_mpi -lboost_serialization
BigRed2CFLAG = -c
BigRed2OFLAG = -o
# General purpose flags.
CC = mpicxx -O3 -g -Wall -fopenmp
LFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_program_options -lboost_mpi -lboost_serialization
CFLAG = -c
OFLAG = -o

all: $(PROG)

install: all
	@echo "Installing $(PROG) into $(homedir) directory"; mv -f $(PROG) $(APPHOME); mkdir $(APPHOME)outfiles

nanoHUB-install:
	. /etc/environ.sh; use -e -r boost-1.62.0-mpich2-1.3-gnu-4.7.2; make CCF=nanoHUB all
	@echo "Installing $(PROG) into $(bindir) directory on Rappture/nanohub"; mv -f $(PROG) $(BIN)

cluster-install:
	module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl; make CCF=BigRed2 all
	@echo "Installing $(PROG) into $(homedir) directory on your local machine"; mv -f $(PROG) $(APPHOME); mkdir $(APPHOME)outfiles

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

main.o: NanoconfinementMd.h
NanoconfinementMd.o: NanoconfinementMd.h utility.h interface.h particle.h vertex.h databin.h control.h functions.h thermostat.h
interface.o: interface.h functions.h
functions.o: functions.h
md.o: particle.h vertex.h interface.h thermostat.h control.h forces.h energies.h functions.h
mdforces.o: forces.h
mdenergies.o: energies.h

clean:
	rm -f *.o $(PROG) 

dataclean:
	rm -f ../outfiles/*.dat ../outfiles/*.xyz ../outfiles/*.lammpstrj ../data/*.dat
	rm -rf $(APPHOME)outfiles

distclean: clean
	rm -f $(BIN)$(PROG) $(APPHOME)$(PROG)

.PHONY: all install clean distclean
