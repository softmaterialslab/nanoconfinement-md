# make cluster-submit-lammps Z=3 p=1 n=-1 c=0.5 d=0.714 S=5000000
# make local-run-lammps Z=3 p=1 n=-1 c=0.5 d=0.714 S=5000000 MPIRUNCMD=mpirun
#This make file builds the sub folder make files
PROG = md_simulation_confined_ions
JOBSCR = nanoconfinement.pbs
TEST = test.pbs
JOBSCRLMP = nanoconfinement-lammps.pbs

BIN = bin
BASE = src
SCRIPT = scripts

Z=3
p=1 
n=-1
c=0.5
d=0.714
S=5000000
NODESIZE=4

MPIRUNCMD = aprun

all:
	@echo "Starting build of the $(BASE) directory";
ifeq ($(CCF),BigRed2)	
	+$(MAKE) -C $(BASE) cluster-install
else ifeq ($(CCF),nanoHUB)
	+$(MAKE) -C $(BASE)
else
	+$(MAKE) -C $(BASE)
endif
	@echo "Ending the build of the $(BASE) directory";
	@echo "installing the $(PROG) into $(BIN) directory"; cp -f $(BASE)/$(PROG) $(BIN)

local-install: create-dirs
	make all	

cluster-install: create-dirs
	make CCF=BigRed2 all

nanoHUB-install: create-dirs
	make CCF=nanoHUB all

create-dirs:
	@echo "Checking and creating needed sub-directories in the $(BIN) directory"
	if ! [ -d $(BIN) ]; then mkdir $(BIN); fi
	if ! [ -d $(BIN)/outfiles ]; then mkdir $(BIN)/outfiles; fi
	if ! [ -d $(BIN)/infiles ]; then mkdir $(BIN)/infiles; fi
	if ! [ -d $(BIN)/data ]; then mkdir $(BIN)/data; fi
	@echo "Directory creation is over."

cluster-submit:
	@echo "Installing jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(JOBSCR) $(BIN)
	+$(MAKE) -C $(BIN) submit

cluster-test-submit:
	@echo "Installing test jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(TEST) $(BIN)
	+$(MAKE) -C $(BIN) test
	
cluster-submit-lammps:
	@echo "Running the preprocessor to create lammps script and input script."
	+$(MAKE) -C $(BIN) run-preprocessor Z=$(Z) p=$(p) n=$(n) c=$(c) d=$(d) S=$(S)
	@echo "Running the preprocessor is over."
	@echo "Installing jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(JOBSCRLMP) $(BIN)
	+$(MAKE) -C $(BIN) submit-lammps
	
run-postprocessor:
	+$(MAKE) -C $(BIN) run-postprocessor Z=$(Z) p=$(p) n=$(n) c=$(c) d=$(d) S=$(S) MPIRUNCMD=$(MPIRUNCMD)
	@echo "Postprocessing is over."

local-run-parallel-lammps:
	@echo "Running the preprocessor to create lammps script and input script."
	+$(MAKE) -C $(BIN) run-preprocessor Z=$(Z) p=$(p) n=$(n) c=$(c) d=$(d) S=$(S) MPIRUNCMD=$(MPIRUNCMD)
	@echo "Running the preprocessor is over."
	+$(MAKE) -C $(BIN) run-local-parallel NODESIZE=$(NODESIZE) MPIRUNCMD=mpirun
	@echo "Lammps simulation is over."
	+$(MAKE) -C $(BIN) run-postprocessor Z=$(Z) p=$(p) n=$(n) c=$(c) d=$(d) S=$(S) MPIRUNCMD=$(MPIRUNCMD)
	@echo "Postprocessing is over."

local-run-lammps:
	@echo "Running the preprocessor to create lammps script and input script."
	+$(MAKE) -C $(BIN) run-preprocessor Z=$(Z) p=$(p) n=$(n) c=$(c) d=$(d) S=$(S) MPIRUNCMD=$(MPIRUNCMD)
	@echo "Running the preprocessor is over."
	+$(MAKE) -C $(BIN) run-local-serial
	@echo "Lammps simulation is over."
	+$(MAKE) -C $(BIN) run-postprocessor Z=$(Z) p=$(p) n=$(n) c=$(c) d=$(d) S=$(S) MPIRUNCMD=$(MPIRUNCMD)
	@echo "Postprocessing is over."


clean: dataclean
	rm -f $(BASE)/*.o
	rm -f $(BASE)/$(PROG)
	rm -f $(BIN)/$(PROG)
	rm -rf $(BIN)/outfiles
	rm -rf $(BIN)/data

dataclean:
	rm -f $(BIN)/outfiles/*.dat $(BIN)/outfiles/*.xyz  $(BIN)/outfiles/*.lammpstrj  $(BIN)/temp/*
	rm -f $(BIN)/data/*.dat $(BIN)/data/*.xyz  $(BIN)/data/*.lammpstrj
	rm -f $(BIN)/*.log
	rm -f $(BIN)/*.pbs

.PHONY: all install cluster-install nanoHUB-install create-dirs cluster-submit cluster-test-submit clean dataclean
