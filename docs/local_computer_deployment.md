---
title: "Deploying Computational Codes on a Local Computer"
keywords: Linux, nanoconfinement, ions, OpenMP/MPI Hybrid
topnav: topnav
hide_sidebar: true
summary: This page has instructions to install the source code for simulating ions in nanoconfinement on a linux computer (here we show specific steps for installation on Ubuntu)
---

## Installing from Source and Building on Linux

The nanoconfinement code is open source and can optionally be built from source and used locally or on computer clusters. The following local build instructions should serve as a reference. 

### Necessary Modules
* Load the necessary modules:
```module load gsl && module load openmpi/3.0.1 && module load boost/1_67_0```
* Make sure to export BOOST_LIBDIR environment variable with location to the lib directory:
```export BOOST_LIBDIR=/opt/boost/gnu/openmpi_ib/lib/```
* Also make sure to export OMP_NUM_THREADS environment variable with maximum threads available in your CPU:
```export OMP_NUM_THREADS=16```

### Install instructions
* Copy or git clone nanoconfinement-md project in to a directory. 
* Go to nanoconfinement-md directory and (cd nanoconfinement-md)
* You should provide the following make command to make the project. This will create the executable and Install the executable (md_simulation_confined_ions) into bin directory (That is nanoconfinement-md/bin)
    * make install
* Next, go to the bin directory:
 ```cd bin```
* Now you are ready to run the executable with aprun command using the following method:
```time mpirun -np 2 -N 16 ./md_simulation_confined_ions -Z 3 -p 1 -n -1 -c 0.5 -d 0.714 -S 1000000```
* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
* Check and compare files (ex: energy.out) inside the ```bin/outfiles``` directory.
* If you want to clean everything and create a new build, use:
```make clean``
* Once the simulation has finished, data and outflies folders will contain the simulation results. You may check final density profile form data folder against the example desity profile provided in nanoconfinement-md/examples folder.
