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
  
void INTERFACE::put_saltions_inside(vector<PARTICLE>& saltion_in, int pz, int nz, double concentration, double diameter, vector<PARTICLE>& ion)
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
  diameter = diameter / unitlength;
  
  // distance of closest approach between the ion and the interface
  double r0_x = 0.5 * lx - 0.5 * diameter;
  double r0_y = 0.5 * ly - 0.5 * diameter;
  double r0_z = 0.5 * lz - 0.5 * diameter;
  
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
    if (x > r0_x - diameter || y > r0_y-diameter || z > r0_z-diameter)		// putting an extra ion diameter length away from interface
      continue;
    bool continuewhile = false;
    for (unsigned int i = 0; i < ion.size() && continuewhile == false; i++)
      if ((posvec - ion[i].posvec).GetMagnitude() <= (0.5*diameter+0.5*ion[i].diameter)) continuewhile = true;
    if (continuewhile == true)
      continue;
    PARTICLE freshion;
    if (saltion_in.size() < total_pions_inside)
      freshion = PARTICLE(int(ion.size())+1,diameter,pz,pz*1.0,1.0,ein,posvec,lx,ly,lz);
    else
      freshion = PARTICLE(int(ion.size())+1,diameter,nz,nz*1.0,1.0,ein,posvec,lx,ly,lz);
    saltion_in.push_back(freshion);		// create a salt ion
    ion.push_back(freshion);			// copy the salt ion to the stack of all ions
  }
  ofstream list_salt_ions_inside("outfiles/salt_ions_inside.xyz", ios::out);
  list_salt_ions_inside << saltion_in.size() << endl;
  list_salt_ions_inside << "salt ions inside" << endl;
  for (unsigned int i = 0; i < saltion_in.size(); i++)
    list_salt_ions_inside << "Si" << setw(15) << saltion_in[i].posvec << endl;
  list_salt_ions_inside.close(); 
  
  gsl_rng_free (rnr);
  
  return;
} 

// discretize interface
void INTERFACE::discretize(double ion_diameter, double f)
{
  // width of the discretization (f is typically 1 or 1/2 or 1/4 or 1/8 ...)
  width = f * ion_diameter;	// in reduced units
  
  // note lx and ly are in units of unitlength; that is they are reduced
  unsigned int nx = int(lx / width);
  unsigned int ny = int(ly / width);
  
  // creating a discretized hard wall interface at z = - l/2
  for (unsigned int j = 0; j < ny; j++)
  {
    for (unsigned int i = 0; i < nx; i++)
    {
      VECTOR3D position = VECTOR3D(-0.5*lx+0.5*width+i*width,-0.5*ly+0.5*width+j*width,-0.5*lz);
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
      VECTOR3D position = VECTOR3D(-0.5*lx+0.5*width+i*width,-0.5*ly+0.5*width+j*width,0.5*lz);
      double area = width * width;
      VECTOR3D normal = VECTOR3D(0,0,1);
      rightplane.push_back(VERTEX(position,area,normal));
    }
  }
  
  ofstream listleftplane("outfiles/leftplane.xyz", ios::out);
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
    listleftplane << k+1 << "  " << "1" << "  " << leftplane[k].posvec << endl;
  listleftplane.close();  
  
  ofstream listrightplane("outfiles/rightplane.xyz", ios::out);
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
    listrightplane << k+1 << "  " << "1" << "  " << rightplane[k].posvec << endl;
  listrightplane.close();  
  
  return;
}