# nanoconfinement-md

Bootstrapping a README file.

Instructions to compile the code on Indiana University BigRed 2 Cluster 
* By default BR2 cluster environment has Cray programming environment module (PrgEnv-cray) loaded 
* Nanoconfinement code has dependency to Boost libraries which requires GNU programming environment
    * Switch modules to GNU - ```module swap PrgEnv-cray PrgEnv-gnu```
* Load latest boost libraries
    * ```module load boost/1.65.0```
* Load GSL libraries
    * ```module load gsl```
* Built the code
    * ```make```
    
What does this code do
* The code enables simulations of ions confined between nanoparticles (NPs) or other material surfaces
    * Length of confinement is of the order of nanometers
* Materials represent nanoparticles (NPs) or biomacromolecules
    * NP surfaces are treated as planar walls due to the large size difference between ions and NPs 
* Users can extract the ionic structure (density profile) for a wide variety of ionic and environmental parameters
* Unpolarized surfaces are assumed and standard molecular dynamics is used to propagrate the dynamics of ions

TODO: 
* Add install instructions on SDSC Comet, TACC Stampede Clusters. 
* Add instructions on how to deploy the code and use it with Jadhao Lab's XSEDE Allocation for nano project.
