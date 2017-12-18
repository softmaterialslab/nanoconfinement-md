// This is main
// This is molecular dynamics (MD) simulation of ions of an electrolyte confined within nanoparticle (NP) surfaces
// NP surfaces are assumed as planar walls
// NP and its environment (typically water) are assumed to have the same dielectric properties
// Problem (output): Compute density profiles of ions trapped within this nanoconfinement

#include "NanoconfinementMd.h"

int main(int argc, char* argv[])
{

  NanoconfinementMd nanoconfinementMd;
  nanoconfinementMd.startSimulation(argc, argv,true);

  return 0;
}
// End of main
