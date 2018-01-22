# This is a makefile.
# Use option -p in CC for profiling with gprof

PROG = md_simulation_confined_ions
BIN = ../bin/
APPHOME = ../
bindir = bin
homedir = home
OBJ = main.o interface.o NanoconfinementMd.o functions.o md.o mdforces.o mdenergies.o

# do not use -fopenmp unless parallelizing properly. 
nanoHUBCC = mpicxx -O3 -g -Wall -fopenmp
nanoHUBLFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_program_options -lboost_mpi -lboost_serialization
nanoHUBCFLAG = -I${BOOST_INCDIR}
nanoHUBOFLAG = -o

BigRed2CC = CC -O3 -g -Wall -fopenmp
BigRed2LFLAG = -lgsl -lgslcblas -lm -lboost_program_options -lboost_mpi -lboost_serialization
BigRed2CFLAG = -c
BigRed2OFLAG = -o

all: $(PROG)

install:
	. /etc/environ.sh; use -e -r boost-1.62.0-mpich2-1.3-gnu-4.7.2; make all
	@echo "Installing $(PROG) into $(bindir) directory on Rappture/nanohub"; mv -f $(PROG) $(BIN)

local-install: all
	@echo "Installing $(PROG) into $(homedir) directory on your local machine"; mv -f $(PROG) $(APPHOME); mkdir $(APPHOME)outfiles

$(PROG) : $(OBJ)
	$(nanoHUBCC) $(nanoHUBOFLAG) $(PROG) $(OBJ) $(nanoHUBLFLAG)

%.o : %.cpp
	$(nanoHUBCC) -c $(nanoHUBCFLAG) $< -o $@

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
	rm -f $(BIN)$(PROG)

.PHONY: all install clean distclean
