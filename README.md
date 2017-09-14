# nanoconfinement-md

Bootstrapping a README file.

Instructions to compile the code on Indiana University BigRed 2 Cluster 
* By default BR2 cluster environment has Cray programming environment module (PrgEnv-cray) loaded. 
* Nanoconfinement code has dependency to Boost libraries which requires GNU programming environment. 
    * Switch modules to GNU - ```module swap PrgEnv-cray PrgEnv-gnu```
* Load latest boost libraries
    * ```module load boost/1.65.0```
* Load GSL libraries
    * ```module load gsl```
* Built the code
    * ```make```

TODO: 
* Add description about what the code does.
* Add install instructions on SDSC Comet, TACC Stampede Clusters. 
* Add instructions on how to deploy the code and use it with Jadhao Lab's XSEDE Allocation for nano project.
