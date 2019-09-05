// This file contains member functions for interface class

#include "interface.h"
#include "functions.h"

void INTERFACE::set_up(double salt_conc_in, double salt_conc_out, int salt_valency_in, int salt_valency_out, double bx, double by, double bz)
{
  // useful combinations of different dielectric constants (inside and outside)
  em = 0.5 * (ein + eout);
  ed = (eout - ein) / (4 * pi);

  // useful length scales signifying competition between electrostatics and entropy
  lB_in = (lB_water * epsilon_water / ein) / unitlength;
  lB_out = (lB_water * epsilon_water / eout) / unitlength;
  if (salt_conc_in != 0)
  {
    inv_kappa_in = (0.257 / (salt_valency_in * sqrt(lB_in * unitlength * salt_conc_in))) / unitlength;
    mean_sep_in = pow(1.2 * salt_conc_in, -1.0/3.0) / unitlength;
  }
  else
  {
    inv_kappa_in = 0;
    mean_sep_in = 0;
  }
  if (salt_conc_out != 0)
  {
    inv_kappa_out = (0.257 / (salt_valency_out * sqrt(lB_out * unitlength * salt_conc_out))) / unitlength;
    mean_sep_out = pow(1.2 * salt_conc_out, -1.0/3.0) / unitlength;
  }
  else
  {
    inv_kappa_out = 0;
    mean_sep_out = 0;
  }

  // simulation box size (in reduced units)
  lx = bx;
  ly = by;
  lz = bz;

  return;
}

void INTERFACE::put_saltions_inside(vector<PARTICLE>& saltion_in, int pz, int nz, double concentration, double positive_diameter_in, double negative_diameter_in, vector<PARTICLE>& ion, unsigned int counterions, int valency_counterion, double counterion_diameter_in, double bigger_ion_diameter)
{
  // establish the number of inside salt ions first
  // Note: salt concentration is the concentration of one kind of ions.
  // also the factor of 0.6 is there in order to be consistent with units.

  double volume_box = lx*ly*lz;

  unsigned int total_pions_inside = int((concentration * 0.6022) * (volume_box * unitlength * unitlength * unitlength));
  if (total_pions_inside % nz !=0)
    total_pions_inside = total_pions_inside - (total_pions_inside % abs (nz));
  unsigned int total_nions_inside = total_pions_inside * pz / abs (nz);
  unsigned int total_saltions_inside = total_nions_inside + total_pions_inside + counterions;

  // express diameter in consistent units
  bigger_ion_diameter = bigger_ion_diameter / unitlength; // the bigger_ion_diameter can be cation or anion depending on their sizes;
  positive_diameter_in = positive_diameter_in / unitlength;
  negative_diameter_in = negative_diameter_in / unitlength;
  counterion_diameter_in = counterion_diameter_in / unitlength;

  // distance of closest approach between the ion and the interface
  // choosing bigger_ion_diameter to define distance of closest approach helps us to avoid overlapping the ions when we generate salt ions inside;
  double r0_x = 0.5 * lx - 0.5 * bigger_ion_diameter;
  double r0_y = 0.5 * ly - 0.5 * bigger_ion_diameter;
  double r0_z = 0.5 * lz - 0.5 * bigger_ion_diameter;

  // UTILITY ugsl;			// utility used for making initial configuration

  const gsl_rng_type * rnT;
  gsl_rng * rnr;

  rnT = gsl_rng_default;
  rnr = gsl_rng_alloc (rnT);

  // unsigned long int s = 23897897;  // gsl_rng_uniform will eventually want a non-negative "long" integer
  gsl_rng_env_setup();

  // if you want to start with the same initial configuration comment three lines below
  // else uncomment them. uncommenting will generate a new randomized distribution of salt ions every run
  // srand((time(0)));                // srand & time are built-in
  // unsigned long int s = random();  // gsl_rng_uniform will eventually
  // gsl_rng_set(rnr,s); // seed the random number generator;

  // generate salt ions inside
  while (saltion_in.size() != total_saltions_inside)
  {
    double x = gsl_rng_uniform(rnr);
    x = (1 - x) * (-r0_x) + x * (r0_x);
    double y = gsl_rng_uniform(rnr);
    y = (1 - y) * (-r0_y) + y * (r0_y);
    double z = gsl_rng_uniform(rnr);
    z = (1 - z) * (-r0_z) + z * (r0_z);
    VECTOR3D posvec = VECTOR3D(x,y,z);
    if (x > r0_x - bigger_ion_diameter || y > r0_y - bigger_ion_diameter  || z > r0_z - bigger_ion_diameter )		// putting an extra ion diameter length away from interface
      continue;
    bool continuewhile = false;
    for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
      if ((posvec - ion[i].posvec).GetMagnitude() <= (0.5*bigger_ion_diameter+0.5*ion[i].diameter)) continuewhile = true;
    if (continuewhile == true)
      continue;
    PARTICLE freshion;
    if (saltion_in.size() < counterions)
      freshion = PARTICLE(int(ion.size())+1,counterion_diameter_in,valency_counterion,valency_counterion*1.0,1.0,ein,posvec,lx,ly,lz);
    else if (saltion_in.size() >= counterions && saltion_in.size() < (total_pions_inside + counterions))
      freshion = PARTICLE(int(ion.size())+1,positive_diameter_in,pz,pz*1.0,1.0,ein,posvec,lx,ly,lz);
    else
      freshion = PARTICLE(int(ion.size())+1,negative_diameter_in,nz,nz*1.0,1.0,ein,posvec,lx,ly,lz);
    saltion_in.push_back(freshion);		// create a salt ion
    ion.push_back(freshion);			// copy the salt ion to the stack of all ions
  }

  mpi::environment env;
  mpi::communicator world;
  if (world.rank() == 0) {
	  string list_salt_ions_insidePath= rootDirectory+"outfiles/salt_ions_inside.xyz";
	  ofstream list_salt_ions_inside(list_salt_ions_insidePath.c_str(), ios::out);
	  list_salt_ions_inside << saltion_in.size() << endl;
	  list_salt_ions_inside << "salt ions inside" << endl;
	  for (unsigned int i = 0; i < saltion_in.size(); i++)
		list_salt_ions_inside << "Si" << setw(15) << saltion_in[i].posvec << endl;
	  list_salt_ions_inside.close();
  }
  gsl_rng_free (rnr);

  return;
}

// discretize interface
void INTERFACE::discretize(double smaller_ion_diameter, double f)
{
  // width of the discretization (f is typically 1 or 1/2 or 1/4 or 1/8 ...)
  width = f * lx;	// in reduced units

  // note lx and ly are in units of unitlength; that is they are reduced
  unsigned int nx = int(lx / width);// + 1;
  unsigned int ny = int(ly / width);// + 1;

  // creating a discretized hard wall interface at z = - l/2
  for (unsigned int j = 0; j < ny; j++)
  {
    for (unsigned int i = 0; i < nx; i++)
    {
      VECTOR3D position = VECTOR3D(-0.5*lx+0.5*smaller_ion_diameter+i*width,-0.5*ly+0.5*smaller_ion_diameter+j*width,-0.5*lz);
      double area = width * width;
      VECTOR3D normal = VECTOR3D(0,0,-1);
      leftplane.push_back(VERTEX(position,area,normal));
    }
  }

  // creating a discretized hard wall interface at z = l/2
  for (unsigned int j = 0; j < ny; j++)
  {
    for (unsigned int i = 0; i < nx; i++)
    {
      VECTOR3D position = VECTOR3D(-0.5*lx+0.5*smaller_ion_diameter+i*width,-0.5*ly+0.5*smaller_ion_diameter+j*width,0.5*lz);
      double area = width * width;
      VECTOR3D normal = VECTOR3D(0,0,1);
      rightplane.push_back(VERTEX(position,area,normal));
    }
  }
  mpi::environment env;
  mpi::communicator world;
  if (world.rank() == 0) {
	  string listleftplanePath= rootDirectory+"outfiles/leftplane.xyz";
	  ofstream listleftplane(listleftplanePath.c_str(), ios::out);
	  listleftplane << "ITEM: TIMESTEP" << endl;
	  listleftplane << 0 << endl;
	  listleftplane << "ITEM: NUMBER OF ATOMS" << endl;
	  listleftplane << leftplane.size() << endl;
	  listleftplane << "ITEM: BOX BOUNDS" << endl;
	  listleftplane << -0.5*lx << "\t" << 0.5*lx << endl;
	  listleftplane << -0.5*ly << "\t" << 0.5*ly << endl;
	  listleftplane << -0.5*lz << "\t" << 0.5*lz << endl;
	  listleftplane << "ITEM: ATOMS index type x y z" << endl;
	  for (unsigned int k = 0; k < leftplane.size(); k++)
		listleftplane << k+1 << "  " << "1" << "  " << leftplane[k].posvec.x <<  "  " <<  leftplane[k].posvec.y <<  "  " << (leftplane[k].posvec.z- 0.5 * smaller_ion_diameter) << endl;
	  listleftplane.close();

	  string listrightplanePath= rootDirectory+"outfiles/rightplane.xyz";
	  ofstream listrightplane(listrightplanePath.c_str(), ios::out);
	  listrightplane << "ITEM: TIMESTEP" << endl;
	  listrightplane << 0 << endl;
	  listrightplane << "ITEM: NUMBER OF ATOMS" << endl;
	  listrightplane << rightplane.size() << endl;
	  listrightplane << "ITEM: BOX BOUNDS" << endl;
	  listrightplane << -0.5*lx << "\t" << 0.5*lx << endl;
	  listrightplane << -0.5*ly << "\t" << 0.5*ly << endl;
	  listrightplane << -0.5*lz << "\t" << 0.5*lz << endl;
	  listrightplane << "ITEM: ATOMS index type x y z" << endl;
	  for (unsigned int k = 0; k < rightplane.size(); k++)
		listrightplane << k+1 << "  " << "1" << "  " << rightplane[k].posvec.x <<  "  " << rightplane[k].posvec.y <<  "  " << (rightplane[k].posvec.z + 0.5 * smaller_ion_diameter) << endl;
	  listrightplane.close();
  }
  return;
}
// creating data file for LAMMPS;
void INTERFACE::generate_lammps_datafile(vector<PARTICLE>& saltion_in, int pz, int nz, vector<PARTICLE>& ion, double smaller_ion_diameter,
  double charge_meshpoint, int counterions, int valency_counterion, double fraction_diameter, double surface_area)
{
  // In Lammps, unit of charge is in reduced unit, where q* = q / (4 pi perm0 sigma epsilon)^1/2;
  string AtomType;
  double ChargeValue;
  double charge_density;
  double unitcharge_lammps = unitcharge / (sqrt(4 * pi * unitlength * pow(10, -9) * 8.854187 * pow(10, -12) * 1.38064852 * pow(10,-23) * room_temperature)); //in reduced unit
  charge_meshpoint = charge_meshpoint * unitcharge_lammps;
  smaller_ion_diameter = smaller_ion_diameter / unitlength;	// should be generalized
  mpi::environment env;
  mpi::communicator world;

  if (world.rank() == 0)
  {
    //we should make sure the total charge of both surfaces and the counter ions are zero;
    if ((unitcharge_lammps * valency_counterion * counterions) + (charge_meshpoint * (leftplane.size() + rightplane.size())) != 0.0)
    {
      charge_meshpoint = -1.0 * (unitcharge_lammps * valency_counterion * counterions) / ((leftplane.size() + rightplane.size()));
      charge_density = ((charge_meshpoint/unitcharge_lammps) * ((unitcharge * pow ((1.0/fraction_diameter), 2.0)))) / surface_area;
      cout  << " In Lammps: charge density of surface is " << charge_density << " Coulomb per squared meter " << endl;
    }

    string InputLammpsPath= rootDirectory+"outfiles/ip.lammps.xyz";
    ofstream listlammps(InputLammpsPath.c_str(), ios::out);
    listlammps << "LAMMPS data file" << endl;
    listlammps << (ion.size() + leftplane.size() + rightplane.size()) << " atoms" << endl;
    listlammps << "3 atom types" << endl; //Type 1 is pz positive charged ions, type 2 is negative charged ions inside the box;
    listlammps << -0.5 * lx << " " << 0.5 * lx<< " " << "xlo xhi" << endl;
    listlammps << -0.5 * ly << " " << 0.5 * ly << " " << "ylo yhi" <<  endl;
    listlammps << (-(0.5 * lz) - pow(10.0,-4))  << " " << ((0.5 * lz) + pow(10.0,-4))  <<  " " << "zlo zhi" << endl;	// if we increase the confinement size by pow(10.0,-4) to make the position of mesh points inside the boundaries;
    listlammps << " " << endl;
    listlammps << "Atoms" << endl;
    listlammps << " " << endl;
    for (unsigned int i = 0; i < ion.size(); i++)
  {
    if (ion[i].valency > 0)
    {
      AtomType = "1";
      ChargeValue = pz * unitcharge_lammps;
    }
    else if (ion[i].valency < 0)
    {
      AtomType = "2";
      ChargeValue = nz * unitcharge_lammps;
    }
    listlammps << i + 1 << "   " << AtomType << "   " << setprecision(10) << ChargeValue << "   " << ion[i].posvec.x << "   " << ion[i].posvec.y << "   " << ion[i].posvec.z << endl;
  }
  for (unsigned int k = 0; k < leftplane.size(); k++)
  {
    listlammps << (k + 1 + ion.size()) << "  " << "3" << "   " << setprecision(10) << charge_meshpoint << "   " << rightplane[k].posvec.x <<  "  " << rightplane[k].posvec.y <<  "  " << rightplane[k].posvec.z << endl;
  }
  for (unsigned int h = 0; h < rightplane.size(); h++)
  {
    listlammps << h + 1 + ion.size() + leftplane.size() << "  " << "3" << "   " << setprecision(10) << charge_meshpoint << "   " << leftplane[h].posvec.x <<  "  " <<  leftplane[h].posvec.y <<  "  " << leftplane[h].posvec.z << endl;
  }
  listlammps.close();
}
  return;
}
