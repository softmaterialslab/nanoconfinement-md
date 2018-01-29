---
title: "Deploying Tools to NanoHUB"
keywords: NanoBio, NCN, NanoHUB
topnav: topnav
hide_sidebar: true
summary: This is a quick cheat sheet based on recent experiences and is by no means a substitute to any of NanoHUB documentations. Use this for quick reference but refer to official documentation for detailed guidance. 
---

NanoHUB facilitates multiple options to enable user interactivity with the tools. The following steps are particularly tuned for NSF funded NCN nodes to deploy tools to NanoHUB Cyber Platform. The steps outlined illustrate one particular approach. These may or may not apply to your use case. Contact NanoHUB support for appropriateness of these instructions to your use case. 

## Set up workspace
* You need to create “affiliated institution” accounts and use your university credentials to sign in to NanoHUB account.
* Tool developers can then submit a ticket to request workspace (Linux desktop in a browser) access
 using the HELP button at the top of any nanoHUB page.

## Convert the project Directory Structure
*  You should convert your project to the above directory structure.
### A typical top-level directory layout
    .
    ├── bin			# Executable will be installed here
    ├── data			# folder for output data
    ├── docs			# Documentation files
    ├── examples		# Example test programs or result files
    ├── middleware		# NanoHub-Rappture specific folder where invoke file is kept
    ├── rappture		# NanoHub-Rappture specific folder where tool.xml and wrapper script are kept	
    ├── scripts			# Script files for cloud/cluster platforms
    ├── src			# Source files of your project
    └── README.md		# Readme file of your project  

## Build Instructions - Makefile
* You need to create a Makefile to bulid and install the project. This Makefile is kept inside src folder. Example Makefile is provided below.
```bash
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
LFLAG = -lgsl -lgslcblas -lm -lboost_program_options -lboost_mpi -lboost_serialization
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

```

## Rappture vs Jupyter 


## Rapptureizing a tool

### GUI Rendering 

### Wrapper Scripts 

### Input Parameters 

### Executing the applications 

### Output Plots

## Deployment Workflow

The following linear workflow is just to get a quick sense. Some of the steps are iterative. 

Status progress: Registered &rarr; Created &rarr; Updated &rarr; Installed &rarr; Approved &rarr; Published

* Register the tool on NanoHUB.org
    * A tool administrator has to manually approve this registration. 
* Providing Source Code
    * The source code is maintained on nanoFORGE svn repository 
    * If the source code is open source and in a public git repository a link to the code can be provided during tool registration. 
    * If the tool is closed source or not using a public git repository, a compressed file of the code can be uploaded during or after the tool registration. 
    * The source code is manually deployed on hub environment by an administrator
* Describe the tool metadata which includes following information 
    * A title for the tool
    * Short description  
    * Abstract or a detailed description about the tool
    * Credits
    * Citations - Typically is auto-generated from NanoHUB link, can use a custom one like a seminal paper about the tool.
    * Funding Sponsors
    * References
    * Publications 
* Updating Source Code
    * Changes to the source code should be communicated to the administrators from the tool status page. Once an update is requested the tool status is changed back to Updated. 
    * NanoHUB tool administrators manually re-deploy the modified code. 
    * The status will be changed to Installed. 
* Once the code is tested and working, the tool owner can approve the tool. 
* The last step is to publish the tool for registered users to see and use it.
* An obvious subset of these will need to be repeated for tool updates. 

