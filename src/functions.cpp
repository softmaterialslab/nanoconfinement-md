// This file contains the routines

#include "functions.h"

// overload out
ostream &operator<<(ostream &os, VECTOR3D vec) {
    os << vec.x << setw(15) << vec.y << setw(15) << vec.z;
    return os;
}

// make bins
void make_bins(vector<DATABIN> &bin, INTERFACE &box, double bin_width) {
    int number_of_bins = int(box.lz / bin_width);
    /*Add two extra bins for contact point densities at both ends*/
    number_of_bins +=2;
    bin.resize(number_of_bins);
    for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
        bin[bin_num].set_up(bin_num, bin_width, box.lx, box.ly, box.lz);
    mpi::environment env;
    mpi::communicator world;
    if (world.rank() == 0) {
        string listbinPath = rootDirectory + "outfiles/listbin.dat";
        ofstream listbin(listbinPath.c_str());
        for (unsigned int num = 0; num < bin.size(); num++)
            listbin << bin[num].n << setw(15) << bin[num].width << setw(15) << bin[num].volume << setw(15)
                    << bin[num].lower << setw(15) << bin[num].higher << endl;
        listbin.close();
    }
    return;
}

// initialize velocities of particles to start simulation
void initialize_particle_velocities(vector<PARTICLE> &ion, vector<THERMOSTAT> &bath) {
    if (bath.size() == 1) {
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].velvec = VECTOR3D(0, 0, 0);                    // initialized velocities
        mpi::environment env;
        mpi::communicator world;
        if (world.rank() == 0)
            cout << "Velocities initialized to 0" << endl;
        return;
    }
    double p_sigma = sqrt(kB * bath[0].T / (2.0 * ion[0].m));        // Maxwell distribution width

    // same random numbers used to generate the gaussian distribution every time. seed is fixed.
    // let me know if you need to change the rnd nums every run.
    UTILITY ugsl;

    for (unsigned int i = 0; i < ion.size(); i++)
        ion[i].velvec = VECTOR3D(gsl_ran_gaussian(ugsl.r, p_sigma), gsl_ran_gaussian(ugsl.r, p_sigma),
                                 gsl_ran_gaussian(ugsl.r, p_sigma));    // initialized velocities
    VECTOR3D average_velocity_vector = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < ion.size(); i++)
        average_velocity_vector = average_velocity_vector + ion[i].velvec;
    average_velocity_vector = average_velocity_vector ^ (1.0 / ion.size());
    for (unsigned int i = 0; i < ion.size(); i++)
        ion[i].velvec = ion[i].velvec - average_velocity_vector;
    return;
}

// make movie
void make_movie(int num, vector<PARTICLE> &ion, INTERFACE &box) {
    mpi::environment env;
    mpi::communicator world;
    if (world.rank() == 0) {
        string outdumpPath = rootDirectory + "outfiles/p.lammpstrj";
        ofstream outdump(outdumpPath.c_str(), ios::app);
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
void compute_n_write_useful_data(int cpmdstep, vector<PARTICLE> &ion, vector<THERMOSTAT> &real_bath, INTERFACE &box,
                                 unsigned int lowerBound,
                                 unsigned int upperBound, vector<double> &ion_energy,
                                 vector<double> &lj_ion_ion, vector<double> &lj_ion_leftdummy,
                                 vector<double> &lj_ion_leftwall, vector<double> &lj_ion_rightdummy,
                                 vector<double> &lj_ion_rightwall) {

    double potential_energy = energy_functional(ion, box, lowerBound, upperBound, ion_energy, lj_ion_ion,
                                                lj_ion_leftdummy, lj_ion_leftwall, lj_ion_rightdummy, lj_ion_rightwall);
    mpi::environment env;
    mpi::communicator world;
    if (world.rank() == 0) {
        string list_temperaturePath = rootDirectory + "outfiles/temperature.dat";
        ofstream list_temperature(list_temperaturePath.c_str(), ios::app);
        string list_energyPath = rootDirectory + "outfiles/energy.dat";
        ofstream list_energy(list_energyPath.c_str(), ios::app);
        list_temperature << cpmdstep << setw(15) << 2 * particle_kinetic_energy(ion) / (real_bath[0].dof * kB)
                         << setw(15) << real_bath[0].T << setw(15) << endl;
        double particle_ke = particle_kinetic_energy(ion);
        double real_bath_ke = bath_kinetic_energy(real_bath);
        double real_bath_pe = bath_potential_energy(real_bath);
        double extenergy = particle_ke + potential_energy + real_bath_ke + real_bath_pe;
        list_energy << cpmdstep << setw(15) << extenergy << setw(15) << particle_ke << setw(15) << potential_energy
                    << setw(15) << particle_ke + potential_energy + real_bath_ke + real_bath_pe << setw(15)
                    << real_bath_ke << setw(15) << real_bath_pe << endl;
        list_temperature.close();
        list_energy.close();
    }
}

// compute density profile of ions
void compute_density_profile(int cpmdstep, double density_profile_samples,
                             vector<double> &meanPositiveionDensity,
                             vector<double> &mean_sq_positiveion_density,
                             vector<double> &meanNegativeionDensity,
                             vector<double> &mean_sq_negativeion_density,
                             vector<PARTICLE> &ion, INTERFACE &box,
                             vector<DATABIN> &bin, CONTROL &cpmdremote) {
    vector<double> sample_positiveion_density;
    vector<double> sample_negativeion_density;

    vector<PARTICLE> positiveion;
    vector<PARTICLE> negativeion;

    for (unsigned int i = 0; i < ion.size(); i++) {
        if (ion[i].valency > 0)
            positiveion.push_back(ion.at(i));
        else if (ion[i].valency < 0)
            negativeion.push_back(ion.at(i));
    }
    mpi::environment env;
    mpi::communicator world;

    bin_ions(positiveion, box, sample_positiveion_density, bin);
    bin_ions(negativeion, box, sample_negativeion_density, bin);

    for (unsigned int b = 0; b < meanPositiveionDensity.size(); b++)
        meanPositiveionDensity.at(b) = meanPositiveionDensity.at(b) + sample_positiveion_density.at(b);
    for (unsigned int b = 0; b < meanNegativeionDensity.size(); b++)
        meanNegativeionDensity.at(b) = meanNegativeionDensity.at(b) + sample_negativeion_density.at(b);
    for (unsigned int b = 0; b < sample_positiveion_density.size(); b++)
        mean_sq_positiveion_density.at(b) =
                mean_sq_positiveion_density.at(b) + sample_positiveion_density.at(b) * sample_positiveion_density.at(b);
    for (unsigned int b = 0; b < sample_negativeion_density.size(); b++)
        mean_sq_negativeion_density.at(b) =
                mean_sq_negativeion_density.at(b) + sample_negativeion_density.at(b) * sample_negativeion_density.at(b);

    // write files
    if ((cpmdstep % cpmdremote.writedensity == 0) && cpmdremote.verbose) {

        std::map<double, std::string> positiveDenistyMap;
        std::map<double, std::string> negativeDensityMap;

        if (world.rank() == 0) {
            char datap[200], datan[200];
            sprintf(datap, "data/_z+_den_%.06d.dat", cpmdstep);
            sprintf(datan, "data/_z-_den_%.06d.dat", cpmdstep);

            string p_density_profile, n_density_profile;
            p_density_profile = rootDirectory + string(datap);
            n_density_profile = rootDirectory + string(datan);

            ofstream outdenp, outdenn;
            outdenp.open(p_density_profile.c_str());
            outdenn.open(n_density_profile.c_str());


            for (unsigned int b = 0; b < meanPositiveionDensity.size(); b++) {
                std::ostringstream stringRow;
                stringRow << bin[b].midPoint * unitlength << setw(15)
                          << meanPositiveionDensity.at(b) / density_profile_samples << endl;

                positiveDenistyMap.insert(
                        std::make_pair(bin[b].midPoint * unitlength, stringRow.str()));

            }

            for (unsigned int b = 0; b < meanNegativeionDensity.size(); b++) {
                std::ostringstream stringRow;
                stringRow << bin[b].midPoint * unitlength << setw(15)
                        << meanNegativeionDensity.at(b) / density_profile_samples << endl;

                negativeDensityMap.insert(
                        std::make_pair(bin[b].midPoint * unitlength, stringRow.str()));

            }

            // Iterate through all elements in std::map to print final denisty plots
            std::map<double, std::string>::iterator itp = positiveDenistyMap.begin();
            while(itp != positiveDenistyMap.end())
            {
                outdenp << itp->second;
                itp++;
            }
           std::map<double, std::string>::iterator itn = negativeDensityMap.begin();
            while(itn != negativeDensityMap.end())
            {
                outdenn << itn->second;
                itn++;
            }

            positiveDenistyMap.clear();
            negativeDensityMap.clear();
            outdenp.close();
            outdenn.close();
        }
    }
    return;
}

void average_errorbars_density(double density_profile_samples, vector<double>& mean_positiveion_density,
                                    vector<double>& mean_sq_positiveion_density,
                                    vector<double>& mean_negativeion_density,
                                    vector<double>& mean_sq_negativeion_density,
                                    vector<PARTICLE>& ion, INTERFACE &box,
                                    vector<DATABIN>& bin, string simulationParams)

{
  // 1. density profile
  vector<double> positiveion_density_profile;
  vector<double> negativeion_density_profile;
  for (unsigned int b = 0; b < mean_positiveion_density.size(); b++)
      positiveion_density_profile.push_back(mean_positiveion_density.at(b) / density_profile_samples);
  for (unsigned int b = 0; b < mean_negativeion_density.size(); b++)
      negativeion_density_profile.push_back(mean_negativeion_density.at(b) / density_profile_samples);

  // 2. error bars
  vector<double> p_error_bar;
  vector<double> n_error_bar;
  for (unsigned int b = 0; b < positiveion_density_profile.size(); b++)
      p_error_bar.push_back(sqrt(1.0 / density_profile_samples) *
                            sqrt(mean_sq_positiveion_density.at(b) / density_profile_samples -
                                 positiveion_density_profile.at(b) * positiveion_density_profile.at(b)));
  for (unsigned int b = 0; b < negativeion_density_profile.size(); b++)
      n_error_bar.push_back(sqrt(1.0 / density_profile_samples) *
                            sqrt(mean_sq_negativeion_density.at(b) / density_profile_samples -
                                 negativeion_density_profile.at(b) * negativeion_density_profile.at(b)));
      // 3. write results
  string p_density_profile, n_density_profile;
  p_density_profile="data/p_density_profile"+simulationParams+".dat";
  n_density_profile="data/n_density_profile"+simulationParams+".dat";
  ofstream list_p_profile(p_density_profile.c_str(), ios::out);
  ofstream list_n_profile(n_density_profile.c_str(), ios::out);
  std::map<double, std::string> positiveDenistyMap;
  std::map<double, std::string> negativeDensityMap;


  for (unsigned int b = 0; b < positiveion_density_profile.size(); b++) {
      std::ostringstream stringRow;
      stringRow << bin[b].midPoint * unitlength << setw(15)<< positiveion_density_profile.at(b) << setw(15) << p_error_bar.at(b) << endl; // change in the z coordinate, counted from leftwall
      positiveDenistyMap.insert(std::make_pair(bin[b].midPoint * unitlength, stringRow.str()));

  }

  for (unsigned int b = 0; b < negativeion_density_profile.size(); b++) {
      std::ostringstream stringRow;
      stringRow << bin[b].midPoint * unitlength << setw(15)<< negativeion_density_profile.at(b) << setw(15) << n_error_bar.at(b)<< endl; // change in the z coordinate, counted from leftwall
      negativeDensityMap.insert(std::make_pair(bin[b].midPoint * unitlength, stringRow.str()));
  }
  // Iterate through all elements in std::map to print final denisty plots
  std::map<double, std::string>::iterator itp = positiveDenistyMap.begin();
  while (itp != positiveDenistyMap.end())
  {
      list_p_profile << itp->second;
      itp++;
  }
  std::map<double, std::string>::iterator itn = negativeDensityMap.begin();
  while (itn != negativeDensityMap.end())
  {
      list_n_profile << itn->second;
      itn++;
  }
  positiveDenistyMap.clear();
  negativeDensityMap.clear();

  list_p_profile.close();
  list_n_profile.close();
return;
}

//Seperate the ljmovie to many data.coords.all.101* files.
//We can skip this function, if in Lammps we write " dump mymovie posneg custom 1000 data.coords.all.101* id  type q  x    y    z"
void output_lammps(vector<PARTICLE>& ion, int &cnt_filename) //cnt_filename shows how many samples are created.
{
  vector<string>lines;
  VECTOR3D posvec;
  string AtomType,ChargeType, Num, line;
  cnt_filename = 0;
  int j = 0;
  int filenumber = 0;
  int header = 9; // There are 9 lines before the atom coordinates start.
  char filename[100];
  ifstream file;
  ofstream outputfile(filename, ios::in);
  file.open("movie_snapshots/ljmovie.xyz");
  while (!file.eof())
  {
    if (j>= 0 && j< header)
    {
      getline(file, line);
      j++;
    }
    if (j >= header && j < (header + ion.size()))
    {
      getline(file,line);
      j++;
      istringstream list(line);
      list >> Num >> AtomType >> ChargeType >> posvec.x >> posvec.y >> posvec.z;
      lines.push_back(line);
    }
    if (j == header + ion.size())
    {
      j=0;
      filenumber = (cnt_filename * 1000);
      sprintf(filename, "movie_snapshots/data.coords.all.%d", filenumber);
      outputfile.open(filename);
      for (vector<string>::iterator it = lines.begin() ; it != lines.end(); ++it)
      {
        outputfile << *it << endl;
      }
      outputfile.close();
      cnt_filename++;
      lines.clear();
    }
  }
  
  return;
}

// Read all data.coords.all.101* files and store them;
//void ReadParticlePositions(vector<PARTICLE>& ion, int i, int data_frequency, int samples, double ion_diameter, double ion_mass, double lx, double ly, double lz, double Charge)
void ReadParticlePositions(vector<PARTICLE>& ion, int cpmdstep, int samples, double diameter, INTERFACE &box)
{
  VECTOR3D posvec;
  string AtomType;
  double Charge;
  int data_frequency = 1000;
  int Num;
  char filename[100];
  int filenumber = cpmdstep * data_frequency;
  sprintf(filename, "movie_snapshots/data.coords.all.%d", filenumber);
  ifstream datafile(filename, ios::in);
  if (!datafile)
  {
    cout << "Position data file " << filenumber << " could not be opened" << endl;
  }
  while (datafile >> Num >> AtomType >> Charge >> posvec.x >> posvec.y >> posvec.z)
  {
    PARTICLE singleljparticle = PARTICLE(int(ion.size())+1,diameter,Charge,0,1.0,box.eout,posvec,box.lx,box.ly,box.lz);
    ion.push_back(singleljparticle);
  }
  return;
}

// compute MD trust factor R
double compute_MD_trust_factor_R(int hiteqm) {
    string inPath = rootDirectory + "outfiles/energy.dat";
    ifstream in(inPath.c_str(), ios::in);
    if (!in) {
        mpi::environment env;
        mpi::communicator world;
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
    mpi::environment env;
    mpi::communicator world;
    if (world.rank() == 0) {
        string outPath = rootDirectory + "outfiles/R.dat";
        ofstream out(outPath.c_str());
        out << "Sample size " << ext.size() << endl;
        out << "Sd: ext, kinetic energy and R" << endl;
        out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;
    }
    return R;
}

// auto correlation function
void auto_correlation_function() {
    string inPath = rootDirectory + "outfiles/for_auto_corr.dat";
    ifstream in(inPath.c_str(), ios::in);
    if (!in) {
        mpi::environment env;
        mpi::communicator world;
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
    mpi::environment env;
    mpi::communicator world;
    if (world.rank() == 0) {
        string outPath = rootDirectory + "outfiles/auto_correlation.dat";
        ofstream out(outPath.c_str());
        for (int i = 0; i < ntau; i++)
            out << i << setw(15) << autocorr[i] / autocorr[0] << endl;

        cout << "Auto correlation function generated" << endl;
    }
    return;
}


// display progress
void ProgressBar(double fraction_completed) {
    int val = (int) (fraction_completed * 100);
    int lpad = (int) (fraction_completed * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% |%.*s%*s|", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
