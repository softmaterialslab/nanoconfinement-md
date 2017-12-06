// Last update: April 5, 2013
// This is main.
// This is MD simulation of an electrolyte confined within planar walls
// There is no dielectric contrast between the walls. This is uniform dielectric medium.
// Problem : Compute density profiles of ions trapped within planar walls

#include "NanoconfinementMd.h"

int main(int argc, char* argv[])
{

  NanoconfinementMd nanoconfinementMd;
  nanoconfinementMd.startSimulation(argc, argv,true);

  return 0;
}
// End of main