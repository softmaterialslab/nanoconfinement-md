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
*  You should convert your project to the following directory structure.
### A typical top-level directory layout
    .
    ├── bin			# Executable will be installed here
    ├── data		# folder for output data
    ├── docs		# Documentation files
    ├── examples		# Example test programs or result files
    ├── middleware		# NanoHub-Rappture specific folder where invoke file is kept
    ├── rappture		# NanoHub-Rappture specific folder where tool.xml and wrapper script are kept	
    ├── scripts		# Script files for cloud/cluster platforms
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
The following workflow explains how to generate a GUI for your application using a XML file and write a wrapper script to pass the input parameters to your program from the rappture GUI and to pass the output from your program to rappture GUI.
### GUI Rendering
* Rappture allows you to render a GUI using a XML file (tool.xml).
* This file should be kept inside the rappture folder in your directory stucture.
* This tool.xml file should have markup tags for input/output parameters, reference to the wrapper script and a description about your application.
* Example tool.xml file is provided below and most markup tags are self explanatory.
```xml
<?xml version="1.0"?>
<run>
    <tool>
        <title>Ions in Nanoconfinement</title>
        <about><!--Describe your project here--></about>
        <command>python @tool/nanoconfinement-wrapper.py @driver</command>
    </tool>
    <input>
        <group id="physical">
            <about>
                <label>Physical</label>
                <description>Physical parameters specific to the material system</description>
            </about>
            <number id="salt_concentration">
                <about>
                    <label>Salt Concentration (c in M)</label>
                    <description>Salt concentration of the solution</description>
                </about>
                <min>0.3</min>
                <max>0.9</max>
                <default>0.5</default>
            </number>
            <integer id="positive_valency">
                <about>
                    <label>Positive Ion Valency (z)</label>
                    <description>Valency (charge) of the positive ions</description>
                </about>
                <min>1</min>
                <max>3</max>
                <default>1</default>
                <preset>
                    <value>1</value>
                    <label>1(monovalent)</label>
                </preset>
                <preset>
                    <value>2</value>
                    <label>2(divalent)</label>
                </preset>
                <preset>
                    <value>3</value>
                    <label>3(trivalent)</label>
                </preset>
            </integer>
            <integer id="negative_valency">
                <about>
                    <label>Negative Ion Valency</label>
                    <description>Valency (charge) of the negative ions</description>
                </about>
                <min>-3</min>
                <max>-1</max>
                <default>-1</default>
                <preset>
                    <value>-1</value>
                    <label>-1</label>
                </preset>
                <preset>
                    <value>-2</value>
                    <label>-2</label>
                </preset>
                <preset>
                    <value>-3</value>
                    <label>-3</label>
                </preset>
            </integer>
            <number id="confinement_length">
                <about>
                    <label>Confinement Length (nm)</label>
                    <description>Length of Nanoscale Confinement</description>
                </about>
                <min>2</min>
                <max>4</max>
                <default>3</default>
            </number>
            <number id="ion_diameter">
                <about>
                    <label>Ion Diameter (nm)</label>
                    <description>Diameter of Ions</description>
                </about>
                <min>0.5</min>
                <max>0.8</max>
                <default>0.714</default>
            </number>
        </group>
        <group id="computing">
            <about>
                <label>Computing</label>
                <description>Computing parameters which effect algorithm</description>
            </about>
            <integer id="simulation_steps">
                <about>
                    <label>Simulation Steps</label>
                    <description>Computing parameters related to the simulation</description>
                </about>
                <min>20000</min>
                <max>5000000</max>
                <default>20000</default>
            </integer>
        </group>
		<image id='ions_distribution'>
			<resize>height=300</resize>
			<current> <!--Base64 version of your image--></current>
		</image>
    </input>
</run>
```
### Wrapper Scripts 
* You can select your preferred programming language to write the wrapper script and here, we explain how to write this wrapper script using programming language Python. 
#### Input Parameters 
* In this wrapper script, you need to get the input parameters from the rappture GUI (generated using tool.xml) using the input element name you have used in the tool.xml. An example is provided below.
```python
io = Rappture.PyXml(sys.argv[1])
salt_concentration = io['input.group(physical).number(salt_concentration).current'].value
```
#### How to link the executable in the wrapper script
* Your executable should be in the bin directory after you do a make-insall. You can link the executable to the wrapper script using a rappture provided library function executeCommand. Follwoing code segment explains how to execute the program. 
```python
try:
     exitStatus,stdOutput,stdError = Rappture.tools.executeCommand(
        ['mpirun','-np','1','md_simulation_confined_ions', '-Z', confinement_length, '-p', positive_valency, '-n', negative_valency, '-c',
         salt_concentration, '-d', ion_diameter, '-S', simulation_steps, '-f', simulation_params, '-v', 'false'], streamOutput=True)
except:
    sys.stderr.write('Error during execution of md_simulation_confined_ions')
    sys.exit(1);
 ```
#### Output Plots

### Executing the applications 

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

