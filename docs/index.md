---
title: "Ions in nanoconfinement"
keywords: nanoconfinement
topnav: topnav
hide_sidebar: true
permalink: index.html
---

## Synopsis

This nanoconfinement code empowers users to simulate ions confined between material surfaces that are nanometers apart, and extract the associated ionic structure. The app facilitates investigations for a wide array of ionic and environmental parameters using standard molecular dynamics (MD) method for unpolarizable surfaces.

## Deploying the tool to NanoHUB
For details on this code has been deployed to NanoHUB Cyber Platform refer [NanoHUB Deployment][nanohub-deployment]

## Using the tool

## Installing from Source

The nanoconfinement code is open source and can optionally be built from source and used locally or on computer clusters. The following cluster build instructions should serve as a reference. 

### Building on IU BigRed 2 Cluster 

* By default BR2 cluster environment has Cray programming environment module (PrgEnv-cray) loaded 
* Nanoconfinement code has dependency to Boost libraries which requires GNU programming environment
    * Switch modules to GNU - ```module swap PrgEnv-cray PrgEnv-gnu```
* Load latest boost libraries
    * ```module load boost/1.65.0```
* Load GSL libraries
    * ```module load gsl```
* Built the code
    * ```make```
