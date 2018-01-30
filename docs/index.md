---
title: "Ions in nanoconfinement"
keywords: nanoconfinement
topnav: topnav
hide_sidebar: true
permalink: index.html
---

## Synopsis

This simulation tool enables researchers to simulate ions confined between nanoparticle (NP) surfaces in aqueous media. Nanoparticles can be synthetic (such as gold NPs) or natural (e.g. proteins) and the length of confinement is of the order of nanometers. Example systems include ion channel proteins of the cell membrane, adsorbed ions near surfaces of porous electrodes, and ions confined by NPs and/or colloidal particles. NP surfaces are assumed to be unpolarizable and are modeled as planar interfaces considering the large size difference between the ions and the NPs. 

The simulation tool facilitates investigations for a wide array of ionic and environmental parameters. Users can extract the ionic structure (density profile) and study its dependence on salt concentration (c), ion valency (z), and other physical attributes. 

You can explore interesting effects by changing the c parameter from 0.3 to 0.9 M. This increase in density leads to crowding of the channel (confinement) with a large number of ions. The effect of symmetry breaking caused by the surfaces is seen: to avoid being pushed by ions from both the sides, an ion prefers the interface over the central region (bulk). The simulations enable the exploration of this effect of ion accumulation near the interface, and make a quantitative assessment of ionic structure in strong confinement.

Another rich avenue to explore is to tune the valency of positive ions (parameter z) from 1 to 3. A positively-charged multivalent ion (+3 Fe or +2 Ca) near an interface is pulled away from the interface by oppositely charged ions with a stronger force relative to the bulk where the symmetry allows for no preferred movement. Thus, stronger electrostatic interactions (as in the case of multivalent ions) tend to cause depletion of the ions from the interface. This code allows you to investigate this depletion effect via accurate computation of the density profiles of ions. 

Effects of changing other physical attributes such as confinement length and ion size are also available to investigate. We invite you to take an inside look at what happens to the self-assembly of ions in these nanoscale channels by investigating the interplay of electrostatic effects and steric (or entropic) effects caused due to confinement, and measuring associated density profiles.

## Deploying the simulation code as Apps on nanoHUB
To deploy the code as a user-friendly app on nanoHUB Cyber Platform, refer to [nanoHUB Deployment](nanohub_deployment)

## Deploying the simulation code on Clusters (BigRed2)
To deploy the code on supercomputing clusters (e.g. BigRed2) to obtain accelerated performance and simulate larger systems for longer times, refer to [Cluster Deployment](cluster_bigred2_deployment)
