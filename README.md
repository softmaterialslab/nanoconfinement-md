# nanoconfinement-md

What does this code do
* The code enables simulations of ions confined between nanoparticles (NPs) or other material surfaces
    * Length of confinement is of the order of nanometers
* Materials represent nanoparticles (NPs) or biomacromolecules
    * NP surfaces are treated as planar walls due to the large size difference between ions and NPs
* Users can extract the ionic structure (density profile) for a wide variety of ionic and environmental parameters
* Unpolarized surfaces are assumed and standard molecular dynamics is used to propagrate the dynamics of ions

Necessary Modules

* Load the necessary modules: module load gsl && module load openmpi/3.0.1 && module load boost/1_67_0
* Make sure to export BOOST_LIBDIR environment variable with location to the lib directory: export BOOST_LIBDIR=/opt/boost/gnu/openmpi_ib/lib/
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
* Note that in current in-house code, you are only able to simulate a system with uncharged surfaces (-i 0.0).

Run the simulation through the LAMMPS:

* In nanoconfinement-md directory, the command to execute the LAMMPS is: make local-run-lammps Z=3 p=1 n=-1 c=0.5 d=0.714 a=0.714 i=0.0 S=1000000 MPIRUNCMD=mpirun LAMMPSEXE=lmp
* This make command creates the input data file (ip.lammps.xyz) through the inhouse code. Then it runs the simulation with LAMMPS. In LAMMPS, you can define charge on the surfaces. For uncharged surfaces: i=0.0. For charged surfaces, you can choose i between zero to -0.01 C/m2. The simulation does not work with charge density more than zero or less than -0.01 C/m-2. 

* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
   * Check and compare files (ex: energy.out) inside the bin/outfiles directory.
   * If you want to clean everything and create a new build, use: ```make clean``
   * Once the simulation has finished, data and outflies folders will contain the simulation results. You may check final density profile form data folder against the example desity profile provided in nanoconfinement-md/examples folder.


For further details please refer to the [documentation](https://softmaterialslab.github.io/nanoconfinement-md/)

## NanoHUB app page:
* https://nanohub.org/tools/nanoconfinement
