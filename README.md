# nanoconfinement-md

What does this code do
* The code enables simulations of ions confined between nanoparticles (NPs) or other material surfaces
    * Length of confinement is of the order of nanometers
* Materials represent nanoparticles (NPs) or biomacromolecules
    * NP surfaces are treated as planar walls due to the large size difference between ions and NPs 
* Users can extract the ionic structure (density profile) for a wide variety of ionic and environmental parameters
* Unpolarized surfaces are assumed and standard molecular dynamics is used to propagrate the dynamics of ions


### Install instructions on local computers

* Copy or git clone nanoconfinement-md project in to a directory.

* Go to nanoconfinement-md/src directory (cd nanoconfinement-md/src/)

* Execute the following command to make the code. This will create and install the executable (md_simulation_confined_ions) into home directory (that is, in nanoconfinement-md/). It will also create necessary output files folder. Some libraries are needed to compile the code (openmpi, boost).

  * make install

* Change directory to the place where the executable is (cd ../). You want to end up in the nanoconfinement-md/ directory.

For further details please refer to the [documentation](https://softmaterialslab.github.io/nanoconfinement-md/) 
