// This file contains the routines 

#include "functions.h"

// overload out
ostream& operator<<(ostream& os, VECTOR3D vec)
{
  os << vec.x << setw(15) << vec.y << setw(15) << vec.z;
  return os;
}

// make bins
void make_bins(vector<DATABIN>& bin, INTERFACE& box, double bin_width)
{
  mpi::environment env;
  mpi::communicator world;

    unsigned int number_of_bins = int(box.lz / bin_width);
    bin.resize(number_of_bins);
    for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
      bin[bin_num].set_up(bin_num, bin_width, box.lx, box.ly, box.lz);
  if (world.rank() == 0) {
    ofstream listbin("outfiles/listbin.dat");
    for (unsigned int num = 0; num < bin.size(); num++)
      listbin << bin[num].n << setw(15) << bin[num].width << setw(15) << bin[num].volume << setw(15) << bin[num].lower
              << setw(15) << bin[num].higher << endl;
    listbin.close();
  }
  return;
}

// initialize velocities of particles to start simulation
void initialize_particle_velocities(vector<PARTICLE>& ion, vector<THERMOSTAT>& bath)
{
  mpi::environment env;
  mpi::communicator world;
  
  if (bath.size() == 1)
  {
    for (unsigned int i = 0; i < ion.size(); i++) 
      ion[i].velvec = VECTOR3D(0,0,0);					// initialized velocities
    if (world.rank() == 0)
	cout << "Velocities initialized to 0" << endl;
    return;
  }
  double p_sigma = sqrt(kB * bath[0].T / (2.0 * ion[0].m));		// Maxwell distribution width
  
  // same random numbers used to generate the gaussian distribution every time. seed is fixed.
  // let me know if you need to change the rnd nums every run.
  UTILITY ugsl;
  
  for (unsigned int i = 0; i < ion.size(); i++) 
    ion[i].velvec = VECTOR3D(gsl_ran_gaussian(ugsl.r,p_sigma), gsl_ran_gaussian(ugsl.r,p_sigma), gsl_ran_gaussian(ugsl.r,p_sigma));	// initialized velocities
  VECTOR3D average_velocity_vector = VECTOR3D(0,0,0);
  for (unsigned int i = 0; i < ion.size(); i++) 
    average_velocity_vector = average_velocity_vector + ion[i].velvec;
  average_velocity_vector = average_velocity_vector ^ (1.0/ion.size());
  for (unsigned int i = 0; i < ion.size(); i++) 
    ion[i].velvec = ion[i].velvec - average_velocity_vector;
  return;
}

// make movie
void make_movie(int num, vector<PARTICLE>& ion, INTERFACE& box)
{
  mpi::environment env;
  mpi::communicator world;
  if (world.rank() == 0) {
    ofstream outdump("outfiles/p.lammpstrj", ios::app);
    outdump << "ITEM: TIMESTEP" << endl;
    outdump << num - 1 << endl;
    outdump << "ITEM: NUMBER OF ATOMS" << endl;
    outdump << ion.size() << endl;
    outdump << "ITEM: BOX BOUNDS" << endl;
    outdump << -0.5 * box.lx << "\t" << 0.5 * box.lx << endl;
    outdump << -0.5 * box.ly << "\t" << 0.5 * box.ly << endl;
    outdump << -0.5 * box.lz << "\t" << 0.5 * box.lz << endl;
    outdump << "ITEM: ATOMS index type x y z" << endl;
    string type;
    for (unsigned int i = 0; i < ion.size(); i++) {
      if (ion[i].valency > 0)
        type = "1";
      else
        type = "-1";
      outdump << setw(6) << i << "\t" << type << "\t" << setw(8) << ion[i].posvec.x << "\t" << setw(8)
              << ion[i].posvec.y << "\t" << setw(8) << ion[i].posvec.z << endl;
    }
    outdump.close();
  }
  return;
}

// compute additional quantities
void compute_n_write_useful_data(int cpmdstep, vector<PARTICLE>& ion, vector<THERMOSTAT>& real_bath, INTERFACE& box)
{
  mpi::environment env;
  mpi::communicator world;
  if (world.rank() == 0) {
  ofstream list_temperature ("outfiles/temperature.dat", ios::app);
  ofstream list_energy ("outfiles/energy.dat", ios::app);
  list_temperature << cpmdstep << setw(15) << 2*particle_kinetic_energy(ion)/(real_bath[0].dof*kB) << setw(15) << real_bath[0].T << setw(15) << endl;
  double particle_ke = particle_kinetic_energy(ion);
  double potential_energy = energy_functional(ion, box);
  double real_bath_ke = bath_kinetic_energy(real_bath);
  double real_bath_pe = bath_potential_energy(real_bath);
  double extenergy = particle_ke + potential_energy + real_bath_ke + real_bath_pe;
  list_energy << cpmdstep << setw(15) << extenergy << setw(15) << particle_ke << setw(15) << potential_energy << setw(15) << particle_ke + potential_energy + real_bath_ke + real_bath_pe << setw(15) << real_bath_ke << setw(15) << real_bath_pe << endl;
  list_temperature.close();
  list_energy.close();

  }
}

// compute MD trust factor R
double compute_MD_trust_factor_R(int hiteqm)
{
  mpi::environment env;
  mpi::communicator world;

    char filename[200];
    sprintf(filename, "outfiles/energy.dat");
    ifstream in(filename, ios::in);
    if (!in) 
    {
	if (world.rank() == 0)
		cout << "File could not be opened" << endl;
      return 0;
    }

    int col1;
    double col2, col3, col4, col5, col6, col7;
    vector<double> ext, ke, pe;
    while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7) {
      ext.push_back(col2);
      ke.push_back(col3);
      pe.push_back(col4);
    }

    double ext_mean = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
      ext_mean += ext[i];
    ext_mean = ext_mean / ext.size();
    double ke_mean = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
      ke_mean += ke[i];
    ke_mean = ke_mean / ke.size();

    double ext_sd = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
      ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
    ext_sd = ext_sd / ext.size();
    ext_sd = sqrt(ext_sd);

    double ke_sd = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
      ke_sd += (ke[i] - ke_mean) * (ke[i] - ke_mean);
    ke_sd = ke_sd / ke.size();
    ke_sd = sqrt(ke_sd);

  double R = ext_sd / ke_sd;
  if (world.rank() == 0) {
    ofstream out("outfiles/R.dat");
    out << "Sample size " << ext.size() << endl;
    out << "Sd: ext, kinetic energy and R" << endl;
    out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;
  }
  return R;
}

// auto correlation function
void auto_correlation_function()
{
  mpi::environment env;
  mpi::communicator world;

    char filename[200];
    sprintf(filename, "outfiles/for_auto_corr.dat");
    ifstream in(filename, ios::in);
    if (!in) 
    {
	if (world.rank() == 0)
		cout << "File could not be opened" << endl;
      return;
    }

    double col1, col2;
    vector<double> n, autocorr;
    while (in >> col1 >> col2)
      n.push_back(col2);

    double avg = 0;
    for (unsigned int j = 0; j < n.size(); j++)
      avg = avg + n[j];
    avg = avg / n.size();

    int ntau = 5000;        // time to which the auto correlation function is computed

    for (int i = 0; i < ntau; i++) {
      double A = 0;
      for (unsigned int j = 0; j < n.size(); j++)
        A = A + n[j + i] * n[j];
      A = A / n.size();
      autocorr.push_back(A - avg * avg);
    }
  if (world.rank() == 0) {
    ofstream out("outfiles/auto_correlation.dat");
    for (int i = 0; i < ntau; i++)
      out << i << setw(15) << autocorr[i] / autocorr[0] << endl;

    cout << "Auto correlation function generated" << endl;
  }
  return;
}
