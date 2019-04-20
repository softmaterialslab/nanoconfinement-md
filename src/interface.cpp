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

void INTERFACE::put_saltions_inside(vector<PARTICLE>& saltion_in, int pz, int nz, double concentration, double positive_diameter_in, double negative_diameter_in, vector<PARTICLE>& ion)
{
  // establish the number of inside salt ions first
  // Note: salt concentration is the concentration of one kind of ions, so for total ions a factor of 2 needs to be multiplied.
  // also the factor of 0.6 is there in order to be consistent with units.

  double volume_box = lx*ly*lz;
  unsigned int total_nions_inside = int((concentration * 0.6022) * (volume_box * unitlength * unitlength * unitlength));
  if (total_nions_inside % pz !=0)
    total_nions_inside = total_nions_inside - (total_nions_inside % pz) + pz;

  unsigned int total_pions_inside = abs(nz) * total_nions_inside / pz;
  unsigned int total_saltions_inside = total_nions_inside + total_pions_inside;

  // express diameter in consistent units
  positive_diameter_in = positive_diameter_in / unitlength;
  negative_diameter_in = negative_diameter_in / unitlength;

  // distance of closest approach between the ion and the interface
  double r0_x = 0.5 * lx - 0.5 * negative_diameter_in;
  double r0_y = 0.5 * ly - 0.5 * negative_diameter_in;
  double r0_z = 0.5 * lz - 0.5 * negative_diameter_in;

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
    if (x > r0_x - negative_diameter_in || y > r0_y-negative_diameter_in || z > r0_z-negative_diameter_in)		// putting an extra ion diameter length away from interface
      continue;
    bool continuewhile = false;
    for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
      if ((posvec - ion[i].posvec).GetMagnitude() <= (0.5*negative_diameter_in+0.5*ion[i].diameter)) continuewhile = true;
    if (continuewhile == true)
      continue;
    PARTICLE freshion;
    if (saltion_in.size() < total_pions_inside)
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
void INTERFACE::discretize(double positive_diameter_in, double f)
{
  // width of the discretization (f is typically 1 or 1/2 or 1/4 or 1/8 ...)
  width = f * lx;	// in reduced units

  // note lx and ly are in units of unitlength; that is they are reduced
  unsigned int nx = int(lx / width) + 1;
  unsigned int ny = int(ly / width) + 1;

  // creating a discretized hard wall interface at z = - l/2
  for (unsigned int j = 0; j < ny; j++)
  {
    for (unsigned int i = 0; i < nx; i++)
    {
      VECTOR3D position = VECTOR3D(-0.5*lx+i*width,-0.5*ly+j*width,-0.5*lz);
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
      VECTOR3D position = VECTOR3D(-0.5*lx+i*width,-0.5*ly+j*width,0.5*lz);
      double area = width * width;
      VECTOR3D normal = VECTOR3D(0,0,1);
      rightplane.push_back(VERTEX(position,area,normal));
    }
  }
  mpi::environment env;
  mpi::communicator world;
  if (world.rank() == 0)
  {
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
		listleftplane << k+1 << "  " << "1" << "  " << leftplane[k].posvec.x <<  "  " <<  leftplane[k].posvec.y <<  "  " << (leftplane[k].posvec.z- 0.5 * positive_diameter_in) << endl;
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
		listrightplane << k+1 << "  " << "1" << "  " << rightplane[k].posvec.x <<  "  " << rightplane[k].posvec.y <<  "  " << (rightplane[k].posvec.z + 0.5 * positive_diameter_in) << endl;
	  listrightplane.close();
  }
  return;
}
// creating data file for LAMMPS;
void INTERFACE::generate_lammps_datafile(vector<PARTICLE>& saltion_in, int pz, int nz, vector<PARTICLE>& ion, double diameter)
{
  // In Lammps, unit of charge is in reduced LJ unit, where q* = q / (4 pi perm0 sigma epsilon)^1/2;
  string AtomType;
  double ChargeValue;
  double Charge = (1.602176634 * pow(10,-19)) / (sqrt(pi * 4.0 * unitlength * pow(10, -9) * 8.854187 * pow(10, -12) * 1.38064852 * pow(10,-23) * room_temperature));
  diameter = diameter / unitlength;
  mpi::environment env;
  mpi::communicator world;

  if (world.rank() == 0)
  {
    string InputLammpsPath= rootDirectory+"outfiles/ip.lammps.xyz";
    ofstream listlammps(InputLammpsPath.c_str(), ios::out);
    listlammps << "LAMMPS data file" << endl;
    listlammps << ion.size() << " atoms" << endl;
    listlammps << "2 atom types" << endl; //Type 1 is pz positive charged ions, type 2 is negative charged ions inside the box;
    listlammps << -0.5 * lx << " " << 0.5 * lx<< " " << "xlo xhi" << endl;
    listlammps << -0.5 * ly << " " << 0.5 * ly << " " << "ylo yhi" <<  endl;
    listlammps << -(0.5 * lz) - (diameter/2.0) << " " << (0.5 * lz) + (diameter/2.0)  <<  " " << "zlo zhi" << endl;
    listlammps << " " << endl;
    listlammps << "Atoms" << endl;
    listlammps << " " << endl;
    for (unsigned int i = 0; i < ion.size(); i++)
    {
      if (ion[i].valency > 0)
      {
        AtomType = "1";
        ChargeValue = pz * Charge;
      }
      else if (ion[i].valency < 0)
      {
        AtomType = "2";
        ChargeValue = nz * Charge;
      }
      listlammps << i + 1 << "   " << AtomType << "   " << ChargeValue << "   " << ion[i].posvec.x << "   " << ion[i].posvec.y << "   " << ion[i].posvec.z << endl;
    }
    listlammps.close();
  }
  return;
}
