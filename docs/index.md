---
title: "Ions in nanoconfinement"
keywords: nanoconfinement
topnav: topnav
hide_sidebar: true
permalink: index.html
---

## Synopsis

This app enables users to simulate ions confined between nanoparticle (NP) surfaces in aqueous media. Nanoparticles can be synthetic (such as gold NPs) or natural (e.g. proteins) and the length of confinement is of the order of nanometers. Example systems include ion channel proteins of the cell membrane, adsorbed ions near surfaces of porous electrodes, and ions confined by NPs and/or colloidal particles. NP surfaces are assumed to be unpolarizable and are modeled as planar interfaces considering the large size difference between the ions and the NPs. 

The app facilitates investigations for a wide array of ionic and environmental parameters. Users can extract the ionic structure (density profile) and study its dependence on salt concentration (c), ion valency (z), and other physical attributes. 

Users can explore interesting effects by changing the c parameter from 0.3 to 0.9 M. This increase in density leads to crowding of the channel (confinement) with a large number of ions. The effect of symmetry breaking caused by the surfaces is seen: to avoid being pushed by ions from both the sides, an ion prefers the interface over the central region (bulk). The app enables users to explore this effect of ion accumulation near the interface, and make a quantitative assessment of ionic structure in strong confinement.

Another rich avenue to explore is to tune the valency of positive ions (parameter z) from 1 to 3. A positively-charged multivalent ion (+3 Fe or +2 Ca) near an interface is pulled away from the interface by oppositely charged ions with a stronger force relative to the bulk where the symmetry allows for no preferred movement. Thus, stronger electrostatic interactions (as in the case of multivalent ions) tend to cause depletion of the ions from the interface. This app empowers users to investigate this depletion effect via accurate computation of the density profiles of ions. 

Effects of changing other physical attributes such as confinement length and ion size are also available for users to explore. We invite users to take an inside look at what happens to the self-assembly of ions in these nanoscale channels by investigating the interplay of electrostatic effects and steric (or entropic) effects caused due to confinement, and measuring associated density profiles.

## Deploying the app to NanoHUB
For details on this code has been deployed to NanoHUB Cyber Platform refer [NanoHUB Deployment](nanohub_deployment)

## Deploying the computational codes to BigRed2

## Installing from Source

The nanoconfinement code is open source and can optionally be built from source and used locally or on computer clusters. The following cluster build instructions should serve as a reference. 

### Building on IU BigRed 2 Cluster 
#### Notes : Following modules are required and automatically loaded in the Makefile
* By default BR2 cluster environment has Cray programming environment module (PrgEnv-cray) loaded 
* Nanoconfinement code has dependency to Boost libraries which requires GNU programming environment
    * Switch modules to GNU - ```module swap PrgEnv-cray PrgEnv-gnu```
* Load latest boost libraries
    * ```module load boost/1.65.0```
* Load GSL libraries
    * ```module load gsl```
#### Install instructions
* Copy or git clone nanoconfinement-md project in to a directory. 
* Go to nanoconfinement-md /src directory and (cd nanoconfinement-md/src/)
* You should provide the following make command to make the project. This will create the executable and Install the executable (md_simulation_confined_ions) into home directory (That is nanoconfinement-md)
    * make cluster-install
* Now you are ready to run the executable with aprun command using any of the following methods : 
##### Using a debug gpu queue.
* Request a debug GPU using the following command. This will provide you 64 processing nodes with 1 hour wall time.
    * qsub -I -l nodes=4:ppn=16,walltime=01:00:00 -q debug_gpu
* change you directory to the place where your executable was :
    * cd nanoconfinement-md
* Load the modules using following command :
    * module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl
* execute the program using following command : 
    * time aprun -n 4 -d 16 ./md_simulation_confined_ions -Z 3 -p 1 -n -1 -c 0.5 -d 0.714 -S 1000000
##### Using a jobscript
* Copy the executable sample script file given in the nanoconfinement-md /scripts folder to the home directory (that is nanoconfinement-md /).
* If you need to change the computational parameter or the physical parameters, you may edit the jobscript file.
* Finally type the following command and it will submit the job.
    * qsub iu_cluster_job_script.pbs 
* You can check the status of the job by typing the following command.
    * qstat -u "IU_username"
* Once the simulation has finished, data and outflies folders will contain the simulation results. You may check final density profile form data folder against the example desity profile provided in nanoconfinement-md /examples/ density_plots folder.

