#This make file builds the sub folder make files
PROG = md_simulation_confined_ions
JOBSCR = iu_cluster_job_script.pbs
TEST = test.pbs

BIN = bin
BASE = src
SCRIPT = scripts

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

install: create-dirs
	make all	

cluster-install: create-dirs
	make CCF=BigRed2 all

nanoHUB-install: create-dirs
	make CCF=nanoHUB all

create-dirs:
	@echo "Checking and creating needed sub-directories in the $(BIN) directory"
	mkdir -p $(BIN)
	mkdir -p $(BIN)/outfiles
	mkdir -p $(BIN)/data
	@echo "Directory creation is over."

cluster-submit:
	@echo "Installing jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(JOBSCR) $(BIN)
	+$(MAKE) -C $(BIN) submit

cluster-test-submit:
	@echo "Installing test jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(TEST) $(BIN)
	+$(MAKE) -C $(BIN) test

clean: dataclean
	rm -f $(BASE)/*.o
	rm -f $(BASE)/$(PROG)
	rm -f $(BIN)/$(PROG)
	rm -rf $(BIN)/outfiles
	rm -rf $(BIN)/data

dataclean:
	rm -f $(BIN)/outfiles/*.dat $(BIN)/outfiles/*.xyz  $(BIN)/outfiles/*.lammpstrj
	rm -f $(BIN)/data/*.dat $(BIN)/data/*.xyz  $(BIN)/data/*.lammpstrj
	rm -f $(BIN)/*.log
	rm -f $(BIN)/*.pbs

.PHONY: all install cluster-install nanoHUB-install create-dirs cluster-submit cluster-test-submit clean dataclean
