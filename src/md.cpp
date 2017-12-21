// This is molecular (particle) dynamics (MD)

#include "particle.h"
#include "vertex.h"
#include "interface.h"
#include "thermostat.h"
#include "control.h"
#include "forces.h"
#include "energies.h"
#include "functions.h"

void
md(vector <PARTICLE> &ion, INTERFACE &box, vector <THERMOSTAT> &real_bath, vector <DATABIN> &bin, CONTROL &mdremote, string &simulationParams) {

    // initialization
    std::vector<vector<VECTOR3D> > forceMatrix(ion.size());
#pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < ion.size(); i++) {
        for (unsigned int j = 0; j < ion.size(); j++)
            forceMatrix[i].push_back(VECTOR3D(0, 0, 0));
    }

    std::vector <VECTOR3D> electrostatic_forces(ion.size(), VECTOR3D(0, 0, 0));
    std::vector <VECTOR3D> lj_ion_ion(ion.size(), VECTOR3D(0, 0, 0));
    std::vector <VECTOR3D> lj_ion_leftdummy(ion.size(), VECTOR3D(0, 0, 0));
    std::vector <VECTOR3D> lj_ion_left_wall(ion.size(), VECTOR3D(0, 0, 0));
    std::vector <VECTOR3D> lj_ion_rightdummy(ion.size(), VECTOR3D(0, 0, 0));
    std::vector <VECTOR3D> lj_ion_right_wall(ion.size(), VECTOR3D(0, 0, 0));



    initialize_particle_velocities(ion, real_bath);    // particle velocities initialized
    for_md_calculate_force(ion, box, 'y',forceMatrix, electrostatic_forces,lj_ion_ion, lj_ion_leftdummy, lj_ion_left_wall, lj_ion_rightdummy, lj_ion_right_wall);        // force on particles and fake degrees initialized
    long double particle_ke = particle_kinetic_energy(ion);// compute initial kinetic energy
    long double potential_energy;
    potential_energy = energy_functional(ion, box);    // compute initial potential energy

    // Output cpmd essentials

    cout << "\n";
    cout << "Propagation of ions using Molecular Dynamics method" << " begins " << endl;
    cout << "Time step in the simulation " << mdremote.timestep << endl;
    cout << "Total number of simulation steps " << mdremote.steps << endl;
    
    if (mdremote.verbose)
    {
      cout << "Initial ion kinetic energy " << particle_ke << endl;
      cout << "Inital potential energy " << potential_energy << endl;
      cout << "Initial system energy " << particle_ke + potential_energy << endl;
      cout << "Chain length (L+1) implementation " << real_bath.size() << endl;
      cout << "Main thermostat temperature " << real_bath[0].T << endl;
      cout << "Main thermostat mass " << real_bath[0].Q << endl;
      cout << "Number of bins used for computing density profiles " << bin.size() << endl;
      cout << "Production begins at " << mdremote.hiteqm << endl;
      cout << "Sampling frequency " << mdremote.freq << endl;
      cout << "Extra computation every " << mdremote.extra_compute << " steps" << endl;
      cout << "Write density profile every " << mdremote.writedensity << endl;
    }

    // for movie
    int moviestart = 1;                    // starting point of the movie

    // for energy
    double energy_samples = 0;

    // for density profile
    vector<double> mean_positiveion_density;            // average density profile
    vector<double> mean_negativeion_density;            // average density profile
    vector<double> mean_sq_positiveion_density;            // average of square of density
    vector<double> mean_sq_negativeion_density;            // average of square of density
    for (unsigned int b = 0; b < bin.size(); b++) {
        mean_positiveion_density.push_back(0.0);
        mean_negativeion_density.push_back(0.0);
        mean_sq_positiveion_density.push_back(0.0);
        mean_sq_negativeion_density.push_back(0.0);
    }
    double density_profile_samples = 0;            // number of samples used to estimate density profile

    long double expfac_real;                // exponential factors useful in velocity Verlet routine

    double percentage=0,percentagePre=-1;

    // Part II : Propagate
    for (int num = 1; num <= mdremote.steps; num++) {
        // INTEGRATOR
        //! begins
        // reverse update of Nose-Hoover chain
        for (int j = real_bath.size() - 1; j > -1; j--)
            update_chain_xi(j, real_bath, mdremote.timestep, particle_ke);
        for (unsigned int j = 0; j < real_bath.size(); j++)
            real_bath[j].update_eta(mdremote.timestep);

        expfac_real = exp(-0.5 * mdremote.timestep * real_bath[0].xi);

        // core loop: velocity --> position --> force --> velocity
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].new_update_velocity(mdremote.timestep, real_bath[0], expfac_real);
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].update_position(mdremote.timestep);
        for_md_calculate_force(ion, box, 'y',forceMatrix, electrostatic_forces,lj_ion_ion, lj_ion_leftdummy, lj_ion_left_wall, lj_ion_rightdummy, lj_ion_right_wall);        // force on particles and fake degrees initialized
        for (unsigned int i = 0; i < ion.size(); i++)
            ion[i].new_update_velocity(mdremote.timestep, real_bath[0], expfac_real);

        // kinetic energies needed to set canonical ensemble
        particle_ke = particle_kinetic_energy(ion);

        // forward update of Nose-Hoover chain
        for (unsigned int j = 0; j < real_bath.size(); j++)
            real_bath[j].update_eta(mdremote.timestep);
        for (unsigned int j = 0; j < real_bath.size(); j++)
            update_chain_xi(j, real_bath, mdremote.timestep, particle_ke);
        //! ends

        // extra computations
        if ((num == 1 || (num % mdremote.extra_compute == 0)) && mdremote.verbose) {
            energy_samples++;
            compute_n_write_useful_data(num, ion, real_bath, box);
            //write_basic_files(num, ion, real_bath, box);	// this is for debug purpose
        }

        // make a movie
        if (num >= moviestart && num % mdremote.moviefreq == 0 && mdremote.verbose)
            make_movie(num, ion, box);

        // compute density profile
        if (num >= mdremote.hiteqm && (num % mdremote.freq == 0)) {
            density_profile_samples++;
            compute_density_profile(num, density_profile_samples, mean_positiveion_density, mean_sq_positiveion_density, mean_negativeion_density, mean_sq_negativeion_density, ion, box, bin, mdremote);
        }

        //percentage calculation
        percentage=roundf(num/(double)mdremote.steps*100 * 10) / 10;
        //percentage output
        if(percentage!=percentagePre){
			if(!mdremote.verbose){
				int progressBarVal=(int) (percentage+0.5);
				printf("=RAPPTURE-PROGRESS=>%d Simulation Running...\n",progressBarVal);
			} else {
                printf("Simulation Completion : %0.1f %%\n",percentage);
            }
			percentagePre=percentage;
        }

    }

    // Part III : Analysis
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
    p_density_profile=rootDirectory+"data/p_density_profile"+simulationParams+".dat";
    n_density_profile=rootDirectory+"data/n_density_profile"+simulationParams+".dat";
    ofstream list_p_profile(p_density_profile.c_str(), ios::out);
    ofstream list_n_profile(n_density_profile.c_str(), ios::out);
    for (unsigned int b = 0; b < positiveion_density_profile.size(); b++)
        list_p_profile << (-0.5 * box.lz + b * bin[b].width) * unitlength << setw(15)
                       << positiveion_density_profile.at(b) << setw(15) << p_error_bar.at(b)
                       << endl; // change in the z coordinate, counted from leftwall
    for (unsigned int b = 0; b < negativeion_density_profile.size(); b++)
        list_n_profile << (-0.5 * box.lz + b * bin[b].width) * unitlength << setw(15)
                       << negativeion_density_profile.at(b) << setw(15) << n_error_bar.at(b)
                       << endl; // change in the z coordinate, counted from leftwall

    string finalConFilePath= rootDirectory+"outfiles/final_configuration.dat";
    ofstream final_configuration(finalConFilePath.c_str());
    for (unsigned int i = 0; i < ion.size(); i++)
        final_configuration << ion[i].posvec << endl;

    if (mdremote.verbose)
    {
      cout << "Number of samples used to compute energy" << setw(10) << energy_samples << endl;
      cout << "Number of samples used to get density profile" << setw(10) << density_profile_samples << endl;
    } 
    cout << "Dynamics of ions simulated for " << mdremote.steps * mdremote.timestep * unittime * 1e9 << " nanoseconds" << endl;

    return;
}
// End of MD routine