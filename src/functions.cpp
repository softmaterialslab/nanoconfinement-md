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
    bin_width = (box.lz / number_of_bins); // To make discretization of bins symmetric, we recalculate the bin_width
    /*Add two extra bins for contact point densities at both ends*/
    number_of_bins += 2;
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
                                 vector<double> &lj_ion_rightwall, vector <double> &coulomb_rightwall, vector <double> &coulomb_leftwall, double charge_meshpoint, int valency_counterion) {

    double potential_energy = energy_functional(ion, box, lowerBound, upperBound, ion_energy, lj_ion_ion,
                                                lj_ion_leftdummy, lj_ion_leftwall, lj_ion_rightdummy, lj_ion_rightwall, coulomb_rightwall, coulomb_leftwall, charge_meshpoint, valency_counterion);
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
            while (itp != positiveDenistyMap.end()) {
                outdenp << itp->second;
                itp++;
            }
            std::map<double, std::string>::iterator itn = negativeDensityMap.begin();
            while (itn != negativeDensityMap.end()) {
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

void average_errorbars_density(double density_profile_samples, vector<double> &mean_positiveion_density,
                               vector<double> &mean_sq_positiveion_density,
                               vector<double> &mean_negativeion_density,
                               vector<double> &mean_sq_negativeion_density,
                               vector<PARTICLE> &ion, INTERFACE &box,
                               vector<DATABIN> &bin, string simulationParams) {
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
    p_density_profile = "data/p_density_profile" + simulationParams + ".dat";
    n_density_profile = "data/n_density_profile" + simulationParams + ".dat";
    ofstream list_p_profile(p_density_profile.c_str(), ios::out);
    ofstream list_n_profile(n_density_profile.c_str(), ios::out);
    std::map<double, std::string> positiveDenistyMap;
    std::map<double, std::string> negativeDensityMap;


    for (unsigned int b = 0; b < positiveion_density_profile.size(); b++) {
        std::ostringstream stringRow;
        stringRow << bin[b].midPoint * unitlength << setw(15) << positiveion_density_profile.at(b) << setw(15)
                  << p_error_bar.at(b) << endl; // change in the z coordinate, counted from leftwall
        positiveDenistyMap.insert(std::make_pair(bin[b].midPoint * unitlength, stringRow.str()));

    }

    for (unsigned int b = 0; b < negativeion_density_profile.size(); b++) {
        std::ostringstream stringRow;
        stringRow << bin[b].midPoint * unitlength << setw(15) << negativeion_density_profile.at(b) << setw(15)
                  << n_error_bar.at(b) << endl; // change in the z coordinate, counted from leftwall
        negativeDensityMap.insert(std::make_pair(bin[b].midPoint * unitlength, stringRow.str()));
    }
    // Iterate through all elements in std::map to print final denisty plots
    std::map<double, std::string>::iterator itp = positiveDenistyMap.begin();
    while (itp != positiveDenistyMap.end()) {
        list_p_profile << itp->second;
        itp++;
    }
    std::map<double, std::string>::iterator itn = negativeDensityMap.begin();
    while (itn != negativeDensityMap.end()) {
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
void output_lammps(vector<PARTICLE> &ion, int &cnt_filename, double data_frequency) //cnt_filename shows how many samples are created.
{
    vector<string> lines;
    VECTOR3D posvec;
    string AtomType, ChargeType, Num, line;
    cnt_filename = 0;
    unsigned int j = 0;
    int filenumber = 0;
    unsigned int header = 9; // There are 9 lines before the atom coordinates start.
    char filename[100];
    ifstream file;

    if (boost::filesystem::exists( "outfiles/electrolyte_movie.xyz" )) {

        ofstream outputfile(filename, ios::in);
        file.open("outfiles/electrolyte_movie.xyz");

        if (boost::filesystem::remove_all("temp") != 0)
            cout << "Pre-existing temp files folder deleted successfully." << endl;

        boost::filesystem::create_directory("temp");

        while (!file.eof()) {
            if (j >= 0 && j < header) {
                getline(file, line);
                j++;
                continue;
            }
            if (j >= header && j < (header + ion.size())) {
                getline(file, line);
                j++;
                istringstream list(line);
                list >> Num >> AtomType >> ChargeType >> posvec.x >> posvec.y >> posvec.z;
                lines.push_back(line);
                continue;
            }
            if (j == header + ion.size()) {
                j = 0;
                filenumber = (cnt_filename * data_frequency);
                sprintf(filename, "temp/data.coords.all.%d", filenumber);
                outputfile.open(filename);
                for (vector<string>::iterator it = lines.begin(); it != lines.end(); ++it) {
                    outputfile << *it << endl;
                }
                outputfile.close();
                cnt_filename++;
                lines.clear();
            }
        }

    } else
        cout << "\noutfiles/electrolyte_movie.xyz not found" << endl;

    return;
}

// Read all data.coords.all.101* files and store them;
//void ReadParticlePositions(vector<PARTICLE>& ion, int i, int data_frequency, int samples, double ion_diameter, double ion_mass, double lx, double ly, double lz, double Charge)
void ReadParticlePositions(vector<PARTICLE> &ion, int cpmdstep, int samples, double diameter, INTERFACE &box, double data_frequency) {
    VECTOR3D posvec;
    string AtomType;
    double Charge;
    int Num;
    char filename[100];
    int filenumber = cpmdstep * data_frequency;
    sprintf(filename, "temp/data.coords.all.%d", filenumber);
    ifstream datafile(filename, ios::in);
    if (!datafile) {
        cout << "Position data file " << filenumber << " could not be opened" << endl;
    }
    while (datafile >> Num >> AtomType >> Charge >> posvec.x >> posvec.y >> posvec.z) {
        PARTICLE singleljparticle = PARTICLE(int(ion.size()) + 1, diameter, Charge, 0, 1.0, box.eout, posvec, box.lx,
                                             box.ly, box.lz);
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

void generateLammpsInputfileForChargedSurface(double ein, int Frequency, int stepsToEqb, int stepsAfterEqb, int extracompute, double timestep, double Positive_diameter_in, double Negative_diameter_in) {

    Positive_diameter_in = Positive_diameter_in / unitlength;
    Negative_diameter_in = Negative_diameter_in / unitlength;
    double Half_positive_diameter_in = 0.5 * Positive_diameter_in;
    double Half_negative_diameter_in = 0.5 * Negative_diameter_in;
    double average_positive_negative_diameter =  0.5 * (Positive_diameter_in + Negative_diameter_in);
    double positive_ion_cutoff = Positive_diameter_in * dcut;
    double negative_ion_cutoff = Negative_diameter_in * dcut;
    double half_positive_ion_cutoff = 0.5 * positive_ion_cutoff;
    double half_negative_ion_cutoff = 0.5 * negative_ion_cutoff;
    double average_negative_positive_cutoff = 0.5 * (negative_ion_cutoff + positive_ion_cutoff);

    /*Replacable variables*/
    string dielectricText = "USERINPUT_DIELECTRIC_CONST";
    string movieFrq = "USERINPUT_MOVIE_FRQ";
    string stepsUpToEQ = "USERINPUT_STEPS_BEFORE_EQ";
    string stepsAfterEQ = "USERINPUT_STEPS_AFTER_EQ";
    string extraComputeSim = "USERINPUT_THERMO_DUMP_FRQ";
    string StepsTime = "USERINPUT_TIMESTEP";

    string positiveIonsDiameter = "USERINPUT_POSITIVE_DIAMETER";
    string negativeIonsDiameter = "USERINPUT_NEGATIVE_DIAMETER";
    string positiveIonsRadius = "POSITIVE_RADIUS";
    string negativeIonsRadius = "NEGATIVE_RADIUS";
    string averagePositiveNegativeIonsDiameter = "AVERAGE_DIAMETER";
    string positiveIonCutoff = "POSITIVE_CUTOFF";
    string negativeIonCutoff = "NEGATIVE_CUTOFF";
    string halfPositiveIonCutoff = "HALF_POS_CUTOFF";
    string halfNegativeIonCutoff = "HALF_NEG_CUTOFF";
    string averageIonCutoff = "CUTOFF_AVERAGE";
    ofstream inputScript("in.lammps", ios::trunc);
    if (inputScript.is_open()) {

        /*Open the template file*/
        string line;
        ifstream inputTemplate("infiles/in.lammps.chargedsurface.template", ios::in);
        if (inputTemplate.is_open()) {
            while (getline(inputTemplate, line)) {
                std::size_t found = line.find(dielectricText);
                if (found != std::string::npos)
                    line.replace(found, dielectricText.length(), std::to_string(ein));

                found = line.find(movieFrq);
                if (found != std::string::npos)
                    line.replace(found, movieFrq.length(), std::to_string(Frequency));

                found = line.find(stepsUpToEQ);
                if (found != std::string::npos)
                    line.replace(found, stepsUpToEQ.length(), std::to_string(stepsToEqb));

                found = line.find(extraComputeSim);
                if (found != std::string::npos)
                    line.replace(found, extraComputeSim.length(), std::to_string(extracompute));

                found = line.find(StepsTime);
                if (found != std::string::npos)
                    line.replace(found, StepsTime.length(), std::to_string(timestep));

                found = line.find(stepsAfterEQ);
                if (found != std::string::npos)
                    line.replace(found, stepsAfterEQ.length(), std::to_string(stepsAfterEqb));

                found = line.find(positiveIonsDiameter);
                if (found != std::string::npos)
                    line.replace(found, positiveIonsDiameter.length(), std::to_string(Positive_diameter_in));

                found = line.find(negativeIonsDiameter);
                if (found != std::string::npos)
                    line.replace(found, negativeIonsDiameter.length(), std::to_string(Negative_diameter_in));

                found = line.find(positiveIonsRadius);
               if (found != std::string::npos)
                   line.replace(found, positiveIonsRadius.length(), std::to_string(Half_positive_diameter_in));

               found = line.find(negativeIonsRadius);
               if (found != std::string::npos)
                   line.replace(found, negativeIonsRadius.length(), std::to_string(Half_negative_diameter_in));

                found = line.find(averagePositiveNegativeIonsDiameter);
                if (found != std::string::npos)
                    line.replace(found, averagePositiveNegativeIonsDiameter.length(), std::to_string(average_positive_negative_diameter));

               found = line.find(positiveIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, positiveIonCutoff.length(), std::to_string(positive_ion_cutoff));

               found = line.find(negativeIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, negativeIonCutoff.length(), std::to_string(negative_ion_cutoff));

               found = line.find(halfPositiveIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, halfPositiveIonCutoff.length(), std::to_string(half_positive_ion_cutoff));

               found = line.find(halfNegativeIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, halfNegativeIonCutoff.length(), std::to_string(half_negative_ion_cutoff));

               found = line.find(averageIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, averageIonCutoff.length(), std::to_string(average_negative_positive_cutoff));

                inputScript << line << endl;
            }
            inputTemplate.close();
        } else cout << "Unable to open the template input script" << endl;
        inputScript.close();
    } else cout << "Unable create a input Script" << endl;

}

void generateLammpsInputfileForUnchargedSurface(double ein, int Frequency, int stepsToEqb, int stepsAfterEqb, int extracompute, double timestep, double Positive_diameter_in, double Negative_diameter_in) {

    Positive_diameter_in = Positive_diameter_in / unitlength;
    Negative_diameter_in = Negative_diameter_in / unitlength;
    double Half_positive_diameter_in = 0.5 * Positive_diameter_in;
    double Half_negative_diameter_in = 0.5 * Negative_diameter_in;
    double average_positive_negative_diameter =  0.5 * (Positive_diameter_in + Negative_diameter_in);
    double positive_ion_cutoff = Positive_diameter_in * dcut;
    double negative_ion_cutoff = Negative_diameter_in * dcut;
    double half_positive_ion_cutoff = 0.5 * positive_ion_cutoff;
    double half_negative_ion_cutoff = 0.5 * negative_ion_cutoff;
    double average_negative_positive_cutoff = 0.5 * (negative_ion_cutoff + positive_ion_cutoff);

    /*Replacable variables*/
    string dielectricText = "USERINPUT_DIELECTRIC_CONST";
    string movieFrq = "USERINPUT_MOVIE_FRQ";
    string stepsUpToEQ = "USERINPUT_STEPS_BEFORE_EQ";
    string stepsAfterEQ = "USERINPUT_STEPS_AFTER_EQ";
    string extraComputeSim = "USERINPUT_THERMO_DUMP_FRQ";
    string StepsTime = "USERINPUT_TIMESTEP";

    string positiveIonsDiameter = "USERINPUT_POSITIVE_DIAMETER";
    string negativeIonsDiameter = "USERINPUT_NEGATIVE_DIAMETER";
    string positiveIonsRadius = "POSITIVE_RADIUS";
    string negativeIonsRadius = "NEGATIVE_RADIUS";
    string averagePositiveNegativeIonsDiameter = "AVERAGE_DIAMETER";
    string positiveIonCutoff = "POSITIVE_CUTOFF";
    string negativeIonCutoff = "NEGATIVE_CUTOFF";
    string halfPositiveIonCutoff = "HALF_POS_CUTOFF";
    string halfNegativeIonCutoff = "HALF_NEG_CUTOFF";
    string averageIonCutoff = "CUTOFF_AVERAGE";
    ofstream inputScript("in.lammps", ios::trunc);
    if (inputScript.is_open()) {

        /*Open the template file*/
        string line;
        ifstream inputTemplate("infiles/in.lammps.unchargedsurface.template", ios::in);
        if (inputTemplate.is_open()) {
            while (getline(inputTemplate, line)) {
                std::size_t found = line.find(dielectricText);
                if (found != std::string::npos)
                    line.replace(found, dielectricText.length(), std::to_string(ein));

                found = line.find(movieFrq);
                if (found != std::string::npos)
                    line.replace(found, movieFrq.length(), std::to_string(Frequency));

                found = line.find(stepsUpToEQ);
                if (found != std::string::npos)
                    line.replace(found, stepsUpToEQ.length(), std::to_string(stepsToEqb));

                found = line.find(extraComputeSim);
                if (found != std::string::npos)
                    line.replace(found, extraComputeSim.length(), std::to_string(extracompute));

                found = line.find(StepsTime);
                if (found != std::string::npos)
                    line.replace(found, StepsTime.length(), std::to_string(timestep));

                found = line.find(stepsAfterEQ);
                if (found != std::string::npos)
                    line.replace(found, stepsAfterEQ.length(), std::to_string(stepsAfterEqb));

                found = line.find(positiveIonsDiameter);
                if (found != std::string::npos)
                    line.replace(found, positiveIonsDiameter.length(), std::to_string(Positive_diameter_in));

                found = line.find(negativeIonsDiameter);
                if (found != std::string::npos)
                    line.replace(found, negativeIonsDiameter.length(), std::to_string(Negative_diameter_in));

                found = line.find(positiveIonsRadius);
               if (found != std::string::npos)
                   line.replace(found, positiveIonsRadius.length(), std::to_string(Half_positive_diameter_in));

               found = line.find(negativeIonsRadius);
               if (found != std::string::npos)
                   line.replace(found, negativeIonsRadius.length(), std::to_string(Half_negative_diameter_in));

                found = line.find(averagePositiveNegativeIonsDiameter);
                if (found != std::string::npos)
                    line.replace(found, averagePositiveNegativeIonsDiameter.length(), std::to_string(average_positive_negative_diameter));

               found = line.find(positiveIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, positiveIonCutoff.length(), std::to_string(positive_ion_cutoff));

               found = line.find(negativeIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, negativeIonCutoff.length(), std::to_string(negative_ion_cutoff));

               found = line.find(halfPositiveIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, halfPositiveIonCutoff.length(), std::to_string(half_positive_ion_cutoff));

               found = line.find(halfNegativeIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, halfNegativeIonCutoff.length(), std::to_string(half_negative_ion_cutoff));

               found = line.find(averageIonCutoff);
               if (found != std::string::npos)
                   line.replace(found, averageIonCutoff.length(), std::to_string(average_negative_positive_cutoff));

                inputScript << line << endl;
            }
            inputTemplate.close();
        } else cout << "Unable to open the template input script" << endl;
        inputScript.close();
    } else cout << "Unable create a input Script" << endl;

}

void NetChargeDensity_ScreenFactor(double charge_density, double bin_width, string simulationParams)
{
  string netcharge_density_profile, screen_factor_profile, errorbars;
  double bin_number, p_density_at_bin_number, n_density_at_bin_number;
  double netcharge_at_bin_number = 0.0;
  double screenfactor_at_bin_number = 0.0;
  char filename[100];
  //open "p_density_profile.dat" and "n_density_profile.dat" files:
  ifstream p_density_file;
  ifstream n_density_file;
  ofstream output_p_density_file(filename, ios::in);
  ofstream output_n_density_file(filename, ios::in);
  p_density_file.open("data/p_density_profile.dat");
  n_density_file.open("data/n_density_profile.dat");
  //create "netcharge_density_profile" and "screen_factor_profile" files:
  netcharge_density_profile = "data/netcharge_density_profile" + simulationParams + ".dat";
  screen_factor_profile = "data/screen_factor_profile" + simulationParams + ".dat";
  ofstream list_netcharge_profile(netcharge_density_profile.c_str(), ios::out);
  ofstream list_screencharge_profile(screen_factor_profile.c_str(), ios::out);

  while ((p_density_file >> bin_number >> p_density_at_bin_number >> errorbars) && (n_density_file >>  bin_number >> n_density_at_bin_number >> errorbars))
  {
    //if p_density_profile.dat and n_density_profile.dat have the same profile, the net charge density will be zero;
    netcharge_at_bin_number = p_density_at_bin_number - n_density_at_bin_number;
    screenfactor_at_bin_number = screenfactor_at_bin_number + (netcharge_at_bin_number * (bin_width * unitlength * pow(10, -9)) * (96485.3399 * 1000 / abs(charge_density)));
    list_netcharge_profile << bin_number << setw(15) << netcharge_at_bin_number << endl;
    // the screen factor is meaningful if there is charge on the walls;
    //otherwis the "list_screencharge_profile" file is empty;
    if (charge_density != 0.0)
    {
      if (bin_number <= 0.0)
      { // we get the screening factor from left wall (-lz/2) to midplane (0);
        list_screencharge_profile << bin_number << setw(15) << screenfactor_at_bin_number << endl;
      }
    }
  }
  list_netcharge_profile.close();
  list_screencharge_profile.close();
  return;
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
