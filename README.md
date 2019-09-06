# nanoconfinement-md

What does this code do
* The code enables simulations of ions confined between nanoparticles (NPs) or other material surfaces
    * Length of confinement is of the order of nanometers
* Materials represent nanoparticles (NPs) or biomacromolecules
    * NP surfaces are treated as planar walls due to the large size difference between ions and NPs
* Users can extract the ionic structure (density profile) for a wide variety of ionic and environmental parameters
* Unpolarized surfaces are assumed and standard molecular dynamics is used to propagrate the dynamics of ions

* THE FOLLOWING INSTRUCTIONS ARE FOR LOCAL INSTALL AND TEST:

Necessary Modules

* Load the necessary modules: module load gsl && module load openmpi/3.0.1 && module load boost/1_67_0
* Also make sure to export OMP_NUM_THREADS environment variable with maximum threads available in your CPU: export OMP_NUM_THREADS=16

Install instructions

* Copy or git clone nanoconfinement-md project in to a directory.
* Go to nanoconfinement-md directory and (cd nanoconfinement-md)
* You should provide the following make command to make the project. This will create the executable and Install the executable (md_simulation_confined_ions) into bin directory (That is nanoconfinement-md/bin)
   * make local-install

* Now, you have two options to run the simulation. You can run the simulation through the in-house code or the LAMMPS. The difference between these two options is the method that we use to calculate the electrostatics energy and force. In in-house code, we calculate the electrostatics energy and force with "charged sheet method". In LAMMPS, the electrostatics energy and force are calculated with Ewald-Summation method.

Run the simulation through the in-house code:

* go to the bin directory: cd bin

* Now you are ready to run the executable with aprun command using the following method: time mpirun -np 2 -N 16 ./md_simulation_confined_ions -Z 3 -p 1 -n -1 -c 0.5 -d 0.714 -a 0.714 -i 0.0 -S 1000000
* In command, you can change many parameters to generate the simulations. Some of these parameters are: -Z length of confinement (nm), -p valency of positive ions, -n valency of negative ions, -c salt concentration (M), -d diameter of positive ion (nm), -a diameter of negative ion (nm), -i charge density on the surface (C.m-2), -S simulation time (steps). The diameter of positive ion can be the same as negative ion (symmetric) or different (asymmetric). You can see some of these examples in nanoconfinement-md/examples. If you want to change more parameters in simulation, please read NanoconfinementMd.cpp file.
* Note that in current in-house code, you are only able to simulate a system with uncharged surfaces (-i 0.0).

Run the simulation through the LAMMPS:

* In nanoconfinement-md directory, the command to execute the LAMMPS is: make local-run-lammps Z=3 p=1 n=-1 c=0.5 d=0.714 a=0.714 i=0.0 S=1000000 MPIRUNCMD=mpirun LAMMPSEXE=lmp
   * Depending on the lammps executable, the following command may work: make local-run-lammps Z=3 p=1 n=-1 c=0.5 d=0.714 a=0.714 i=0.0 S=1000000 MPIRUNCMD=mpirun LAMMPSEXE=lmp_g++
* This make command creates the input data file (ip.lammps.xyz) through the inhouse code. Then it runs the simulation with LAMMPS. In LAMMPS, you can define charge on the surfaces. For uncharged surfaces: i=0.0. For charged surfaces, you can choose i between zero to -0.01 C/m2. The simulation does not work with charge density more than zero or less than -0.01 C/m-2.

* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
   * Check and compare files inside the bin/outfiles and bin/data directories. In bin/data directory, you can compare the averaged density profiles of cation and anion (with error bars). If you ran the simulation with in-house code, the energy values are stored in outfiles/energy.out (1st column: step number, 3rd: kinetic energy, 4th: potential energy). If you ran the simulation with LAMMPS code, the energy values are stored in outfiles/thermo.dat (1st column: step number, 6th: kinetic_energy per particle, 7th: potential_energy per particle). If you wish to compare the results of energy profiles from inhouse code (energy.dat) and LAMMPS (thermo.dat), in thermo.dat you should multiply the energy values by total number of atoms (the total number of atoms is printed in 2nd line in ip.lammps.xyz).

   * If you want to clean everything and create a new build, use: ```make clean```
   * Once the simulation has finished, data and outflies folders will contain the simulation results. You may check final density profile form data folder against the example desity profile provided in nanoconfinement-md/examples folder.


For further details please refer to the [documentation](https://softmaterialslab.github.io/nanoconfinement-md/)

## NanoHUB app page:
* https://nanohub.org/tools/nanoconfinement
