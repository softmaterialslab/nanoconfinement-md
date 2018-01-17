# This is a makefile.
# Use option -p in CC for profiling with gprof

PROG = md_simulation_confined_ions
BIN = ../bin/
HOME = ../
bindir = bin
homedir = home
OBJ = main.o interface.o NanoconfinementMd.o functions.o md.o mdforces.o mdenergies.o

# do not use -fopenmp unless parallelizing properly. 
CC = mpich -O3 -g -Wall -fopenmp

LFLAG = -lgsl -lgslcblas -lm -lboost_program_options -lboost_mpi -lboost_serialization

CFLAG = -c

OFLAG = -o

all: $(PROG)

install: all
	@echo "Installing $(PROG) into $(bindir) directory on Rappture/nanohub"; mv -f $(PROG) $(BIN)

local-install: all
	@echo "Installing $(PROG) into $(homedir) directory on your local machine"; mv -f $(PROG) $(HOME); mkdir $(HOME)outfiles

$(PROG) : $(OBJ)
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LIBS) $(LFLAG)

%.o : %.cpp
	$(CC) -c $(LFLAG) $< -o $@

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
	rm -f ../outfiles/*.dat ../outfiles/*.xyz ../outfiles/*.lammpstrj ../data/*.dat $(BIN)$(PROG) $(HOME)$(PROG); rmdir $(HOME)outfiles

distclean: clean dataclean

.PHONY: all install clean distclean
