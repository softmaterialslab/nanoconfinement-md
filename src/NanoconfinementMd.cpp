// This is main.
// This is MD simulation of an electrolyte confined within planar walls
// Uniform dielectric media are assumed
// Problem : Compute density profiles of ions trapped within planar walls
/* Useful studies :
		     1. Role of valency of ions
		     2. Role of varying salt concentration
		     3. Useful parameters via boost: -Z 3 -p 1 -n -1 -c 0.5 -d 0.714 -S 10000
*/
#include "NanoconfinementMd.h"

double unitlength;
double unittime;
double scalefactor;

int NanoconfinementMd::startSimulation(int argc, char *argv[], bool paraMap) {
    // Electrostatic system variables
    double bx, by, bz;        // lengths of the box
    double ein;            // permittivity of inside medium
    double eout;            // permittivity of outside medium
    int pz_in;            // positive valency of ions inside
    int nz_in;            // negative valency of ions inside
    double salt_conc_in;        // salt concentration outside	(enter in M)
    double positive_diameter_in;    // positive ion diameter
    double negative_diameter_in;   // negative ion diameter
    double counterion_diameter_in = 0.0; // counterion ion diameter; (counterion ions assumed to be positive diameter and surfaces are negatively charged)
    double T;            // temperature at which the system of ions is
    double charge_meshpoint = 0.0; // charge on mesh points to create uniform charge density on surface
    double charge_density; // charge density on surface
    int valency_counterion;
    unsigned int counterions = 0.0;         //number of counter ions
    double total_surface_charge = 0.0; //total charge on the surface (in unit of electron charge)
    double surface_area = 0.0; // area of surface
    double number_meshpoints = 0.0; // number of mesh points on the surface
    double smaller_ion_diameter = 0.0; // the salt ion with smaller size;
    double bigger_ion_diameter = 0.0; // the salt ion with bigger size;


    // Simulation related variables
    double fraction_diameter;        // fraction that multiplies the diameter to generate the discretization width for the interface
    double Q;                // thermostat mass required to generate canonical ensemble
    unsigned int chain_length_real;    // Nose Hoover thermostat chain length for particles
    double bin_width;            // width of the bins used to compute density profiles
    bool lammps = true; //if it is false, do the simulation with c++ code, if it is true do it with lammps;
    bool lammpsPreprocessing = true;
    bool screen = false; //if it is false, do not calculate screen factor
    CONTROL mdremote;            // remote control for md
    string config_file = "input_config.cfg";         //Configuration file path holder
    string simulationParams;

    // Different parts of the system
    vector<PARTICLE> saltion_in;        // salt ions inside
    vector<PARTICLE> ion;            // all ions in the system

    // Analysis
    vector<DATABIN> bin;            // bins
    mpi::environment env;
    mpi::communicator world;

    // Get input values from the user
    //options_description desc("Usage:\nrandom_mesh <options>");
    // Declare a group of options that will be
    // allowed only on command line
    options_description generic("Generic options");
    generic.add_options()
            ("help,h", "print usage message")
            ("bx,X", value<double>(&bx)->default_value(15.3153),
             "box length in x direction in nanometers")    // enter in nanometers
            ("by,Y", value<double>(&by)->default_value(15.3153),
             "box length in y direction in nanometers")    // enter in nanometers
            ("epsilon_in,e", value<double>(&ein)->default_value(80),
             "dielectric const inside")        // must have ein = eout
            ("epsilon_out,E", value<double>(&eout)->default_value(80),
             "dielectric const outside")        // must have ein = eout
            ("fraction_diameter,g", value<double>(&fraction_diameter)->default_value(1/28.0),
             "for interface discretization width")    // examples: 1/28 and 1/50
            ("thermostat_mass,Q", value<double>(&Q)->default_value(1.0), "thermostat mass")
            ("chain_length_real,L", value<unsigned int>(&chain_length_real)->default_value(5),
             "chain length for real system: enter L+1 if you want L thermostats")
            ("bin_width,B", value<double>(&bin_width)->default_value(0.05), "bin width (reduced units)")// in reduced units
            ("md_timestep,T", value<double>(&mdremote.timestep)->default_value(0.001), "time step used in md (reduced units)")
            ("md_eqm,P", value<int>(&mdremote.hiteqm)->default_value(100000), "production begin (md)")
            ("md_freq,F", value<int>(&mdremote.freq)->default_value(100), "sample frequency (md)")
            ("md_extra_compute,x", value<int>(&mdremote.extra_compute)->default_value(10000),
             "compute additional (md)")
            ("md_writedensity,w", value<int>(&mdremote.writedensity)->default_value(100000), "write density files")
            ("md_movie_freq,m", value<int>(&mdremote.moviefreq)->default_value(10000), "compute additional (md)")
            ("simulation_params,f", value<string>(&simulationParams)->default_value(""), "Simulation parameters")
            ("verbose,v", value<bool>(&mdremote.verbose)->default_value(true), "verbose true: provides detailed output");
            //("config,conf", value<string>(&config_file)->default_value("input_config.cfg"),
            //"name of a file of a configuration.")

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    //-Z 3 -p 1 -n -1 -c 0.5 -d 0.714 -S 5000000
    options_description config("Configuration");
    config.add_options()
            ("confinement_length,Z", value<double>(&bz)->default_value(3.0),
             "box length in z direction in nanometers")        // enter in nanometers
            ("positive_valency,p", value<int>(&pz_in)->default_value(1), "positive valency inside")
            ("negative_valency,n", value<int>(&nz_in)->default_value(-1), "negative valency inside")
            ("salt_concentration,c", value<double>(&salt_conc_in)->default_value(0.50), "salt concentration inside (M)")
            ("positive_diameter,d", value<double>(&positive_diameter_in)->default_value(0.474),
             "positive ion diameter inside (nm)")        // enter in nanometers
             ("negative_diameter,a", value<double>(&negative_diameter_in)->default_value(0.627),
              "negative ion diameter inside (nm)")        // enter in nanometers
             ("chargedensity_surface,i", value<double>(&charge_density)->default_value(0.0),
               "charge density on surface (C/m2)")        // enter in c/m2
            ("simulation_steps,S", value<int>(&mdremote.steps)->default_value(5000000), "steps used in md")
            ("lammps,J", value<bool>(&lammps)->default_value(false), "LAMMPS (true LAMMPS; false MD)")
            ("lammpsPreprocessing,j", value<bool>(&lammpsPreprocessing)->default_value(true), "LAMMPS Preprocessing/Postprocessing (true Preprocessing; false Postprocessing)");

    options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    variables_map vm;
    store(parse_command_line(argc, argv, cmdline_options), vm);
    notify(vm);

    if (world.rank() == 0 && (!lammps))
    {
      cout << "\nSimulation begins\n";
    }

    if (world.rank() == 0) {
        if (mdremote.verbose) {
            cout << "For help with the menu, type ./md_simulation_confined_ions -h" << endl;
            cout << "The default app simulates a total of 424 ions" << endl;
        }

        cout << "MPI/OpenMP hybrid acceleration is built in this app" << endl;
        cout << "-----------------------------------------------------" << endl;
    }

    ifstream ifs(config_file.c_str());
    if (!ifs && mdremote.verbose && world.rank() == 0) {
        cout << "can not open config file: " << config_file << "\n";
        //return 0;
    } else {
        store(parse_config_file(ifs, config), vm);
        notify(vm);
    }
    if (world.rank() == 0) {
        if (vm.count("help")) {
            std::cout << cmdline_options << "\n";
            return 0;
        }
    }

    //X and Y mapping
    if (paraMap) {
	if (positive_diameter_in <= negative_diameter_in)
	{
	    unitlength = positive_diameter_in;
	    smaller_ion_diameter = positive_diameter_in;
	    bigger_ion_diameter = negative_diameter_in;
	}
	else
	{
	    unitlength = negative_diameter_in;
	    smaller_ion_diameter = negative_diameter_in;
	    bigger_ion_diameter = positive_diameter_in;
	}

        unittime = sqrt(unitmass * unitlength * pow(10.0, -7) * unitlength / unitenergy);
        scalefactor = epsilon_water * lB_water / unitlength;
        bx = sqrt(212 / 0.6022 / salt_conc_in / bz);
        //bx = 10.0;
        by = bx;

        if ( charge_density < -0.01 || charge_density > 0.0) //we can choose charge density on surface between 0.0 (uncharged surfaces)  to -0.01 C/m2.
        {
          cout << "\ncharge density on the surface must be between zero to -0.01 C/m-2; aborting\n";
          return 0;
        }


        valency_counterion = 1; //pz_in;
        counterion_diameter_in = positive_diameter_in;
        surface_area = bx * by * pow(10.0,-18);// in unit of squared meter;
        number_meshpoints =  pow ((1.0/fraction_diameter), 2.0);
        charge_meshpoint = (charge_density * surface_area) / (unitcharge * number_meshpoints); //in unit of electron charge;
        total_surface_charge = charge_meshpoint * number_meshpoints; //in unit of electron charge;
        counterions =  2.0 * (int( abs (total_surface_charge)/valency_counterion)); // there are two charged surfaces, we multiply the counter ions by two;

        //we should make sure the total charge of both surfaces and the counter ions are zero;
        if (((valency_counterion * counterions) + (total_surface_charge * 2.0 )) != 0)
        { //we distribute the extra charge to the mesh points to make the system electroneutral; then we recalculate the charge density on surface;
          charge_meshpoint = -1.0 * (valency_counterion * (counterions/2.0)) / (number_meshpoints);
          total_surface_charge = -1.0 * (valency_counterion * (counterions/2.0)); //we recalculate the total charge on teh surface;
          charge_density = (total_surface_charge * unitcharge) / surface_area; //in unit of Coulomb per squared meter;
        }

        if (mdremote.steps < 100000) {      // minimum mdremote.steps is 20000
            mdremote.hiteqm = (int)(mdremote.steps*0.1);
            mdremote.writedensity =(int)(mdremote.steps*0.1);
            mdremote.extra_compute = (int)(mdremote.steps*0.01);
            mdremote.moviefreq = (int)(mdremote.steps*0.001);
        }
        else {
            mdremote.hiteqm = (int) (mdremote.steps * 0.2);
            mdremote.writedensity = (int) (mdremote.steps * 0.1);
            mdremote.extra_compute = (int) (mdremote.steps * 0.01);
            mdremote.moviefreq = (int) (mdremote.steps * 0.001);
        }
    }

    // Set up the system
    T = 1;        // set temperature (in reduced units; see utility.h)
    INTERFACE box = INTERFACE(VECTOR3D(0, 0, 0), ein, eout); // interface, z planes hard walls; rest periodic boundaries
    box.set_up(salt_conc_in, 0, pz_in, 0, bx / unitlength, by / unitlength, bz / unitlength);
    box.put_saltions_inside(saltion_in, pz_in, nz_in, salt_conc_in, positive_diameter_in, negative_diameter_in, ion, counterions, valency_counterion, counterion_diameter_in, bigger_ion_diameter);
    make_bins(bin, box, bin_width);    // set up bins to be used for computing density profiles
    /*This is to get contact point densities*/
    double leftContact = -0.5 * box.lz + 0.5 * ion[0].diameter - 0.5 * bin[0].width;
    double rightContact = 0.5 * box.lz - 0.5 * ion[0].diameter - 0.5 * bin[0].width;
    bin[bin.size() - 1].lower = leftContact;
    bin[bin.size() - 2].lower = rightContact;
    bin[bin.size() - 1].higher = leftContact + bin[0].width;
    bin[bin.size() - 2].higher = rightContact + bin[0].width;
    bin[bin.size() - 1].midPoint = 0.5 * (bin[bin.size() - 1].lower + bin[bin.size() - 1].higher);
    bin[bin.size() - 2].midPoint = 0.5 * (bin[bin.size() - 2].lower + bin[bin.size() - 2].higher);

    vector<double> initial_density;
    bin_ions(ion, box, initial_density, bin);    // bin the ions to get initial density profile

    box.discretize(smaller_ion_diameter / unitlength, fraction_diameter, charge_meshpoint);

    if (world.rank() == 0) {
        // output to screen the parameters of the problem
        cout << "Ions are confined by nanomaterial surfaces in aqueous solvent" << endl;
        cout << "Material surfaces modeled as thin planar interfaces" << endl;
        cout << "Solvent modeled as implicit media" << endl;
        cout << "Ions modeled as finite-size, soft spheres" << endl;
        cout << "Dielectric constant of water " << epsilon_water << endl;
        cout << "Unit of length is " << unitlength << " nanometers"
             << endl; // half Bjerrum length; close to Na ion radius
        cout << "Unit of mass is " << unitmass << " grams" << endl; // mass of sodium atom in CGS; in grams
        cout << "Unit of energy is " << unitenergy << " ergs (CGS)" << endl; // kB room_T
        cout << "Unit of time is " << unittime << " seconds"
             << endl; // reduced units (LJ); you can derive from above 3 units.
        cout << "Simulation employs reduced units to measure physical quantities" << endl;
        cout << "Simulation box dimensions (in reduced units) x | y | z " << "  " << box.lx << " | " << box.ly << " | "
             << box.lz << endl;
        cout << "Box dimensions (in nanometers) x | y | z " << "  " << bx << " | " << by << " | " << bz << endl;
        cout << "Note: periodic boundaries in x and y direcitons " << endl;
        cout << "Permittivity inside the confinement (channel) " << box.ein << endl;
        cout << "Permittivity outside " << box.eout << endl;
        cout << "Dielectric contrast across interfaces " << 2 * (box.eout - box.ein) / (box.eout + box.ein) << endl;
        cout << "Positive (+) ion valency " << pz_in << endl;
        cout << "Negative (-) ion valency " << nz_in << endl;
        cout << "Valency of counter ions is " << valency_counterion << endl;
        cout << "positive ion diameter " << positive_diameter_in / unitlength << endl;
        cout << "negative ion diameter " << negative_diameter_in / unitlength << endl;
        cout << "counter ion diameter " << counterion_diameter_in / unitlength << endl;
        cout << "In MD, charge density " << charge_density << " Coulomb per squared meter" << endl;
        cout << "Ion (salt) concentration (c) inside " << salt_conc_in << " M" << endl;
        cout << "Note: we define c = total number of positive ions / volume" << endl;
        cout << "Debye length " << box.inv_kappa_in << endl;
        cout << "Mean separation between ions " << box.mean_sep_in << endl;
        cout << "Temperature (in Kelvin) " << room_temperature << endl;

        if (mdremote.verbose) {
            cout << "Binning width (uniform) " << bin[0].width << endl;
            cout << "Number of bins " << bin.size() << endl;
            cout << "Number of points discretizing the left and right planar walls/interfaces/surfaces "
                 << box.leftplane.size() << "  " << box.rightplane.size() << endl;
        }

        // write to files

        // initial density
        string density_profilePath = rootDirectory + "outfiles/initial_density_profile.dat";
        ofstream density_profile(density_profilePath.c_str(), ios::out);
        for (unsigned int b = 0; b < initial_density.size(); b++)
            density_profile << bin[b].midPoint << setw(15) << initial_density.at(b) << endl;
        density_profile.close();

        // check point
        double totalions = 0;
        for (unsigned int b = 0; b < initial_density.size(); b++) {
            if (bin[b].lower == leftContact || bin[b].higher == rightContact)
                continue;
            totalions += initial_density.at(b) * bin[b].volume;
        }
        int totalpions = 0, totalnions = 0;
        for (unsigned int i = 0; i < ion.size(); i++) {
            if (ion[i].valency > 0)
                totalpions += 1;
            else if (ion[i].valency < 0)
                totalnions += 1;
        }

        if (mdremote.verbose) {
            cout << "Number of ions " << totalions << endl;
            cout << "Number of positive ions " << totalpions << endl;
            cout << "Number of negative ions " << totalnions << endl;
            cout << "Number of counter ions " << counterions << endl;
        }

        if (box.total_charge_inside(ion) + (total_surface_charge * 2.0 ) == 0)
            cout << "System simulated is electroneutral-- total charge inside is 0" << endl;
        else {
            cout << "System not electroneutral; aborting" << endl;
            cout << "Total charge inside the confinement " << box.total_charge_inside(ion) - (total_surface_charge * 2.0) << endl;
            return 0;
        }

        int numOfNodes = world.size();

#pragma omp parallel default(shared)
        {
            if (omp_get_thread_num() == 0) {
                printf("The app comes with MPI and OpenMP (Hybrid) parallelization)\n");
                printf("Number of MPI processes used %d\n", numOfNodes);
                printf("Number of OpenMP threads per MPI process %d\n", omp_get_num_threads());
                printf("Make sure that number of grid points and ions is greater than %d\n",
                       omp_get_num_threads() * numOfNodes);
            }
        }
    }

    // prepare for md : make real baths
    vector<THERMOSTAT> real_bath;
    if (chain_length_real == 1)
        real_bath.push_back((THERMOSTAT(0, T, 3 * ion.size(), 0.0, 0, 0)));
    else {
        real_bath.push_back((THERMOSTAT(Q, T, 3 * ion.size(), 0, 0, 0)));
        while (real_bath.size() != chain_length_real - 1)
            real_bath.push_back((THERMOSTAT(Q / (3 * ion.size()), T, 1, 0, 0, 0)));
        real_bath.push_back((THERMOSTAT(0, T, 3 * ion.size(), 0.0, 0, 0)));
        // final bath is dummy bath (dummy bath always has zero mass)
    }
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

    // Simulation using Molecular Dynamics
    if (!lammps) {

        md(ion, box, real_bath, bin, mdremote, simulationParams, charge_meshpoint, valency_counterion, screen);

        // Post simulation analysis (useful for short runs, but performed otherwise too)
        if (world.rank() == 0) {
            if (mdremote.verbose) {
                cout << "MD trust factor R (should be < 0.05) is " << compute_MD_trust_factor_R(mdremote.hiteqm) << endl;
                // Perform the following calculation when testing for how frequently you should sample data to ensure the samples are decorrelated
                //auto_correlation_function();
            }

            cout << "Converged results for ionic densities expected with ~ 1500 nanoseconds of ion dynamics" << endl;
            cout << "Accordingly, we recommend using ~ 1000,000 simulation steps to obtain smoother, converged profiles "
                 << endl;
            cout << "Simulation ends \n\n";
        }

    } else {

        if (world.rank() == 0) {

            if (lammpsPreprocessing) {

                cout << "Lammps Preprocessing started." << endl;
                if (charge_meshpoint == 0)
                {
                  box.generate_lammps_datafile_unchargedsurface(saltion_in, pz_in, nz_in, ion);
                  generateLammpsInputfileForUnchargedSurface(ein, mdremote.freq, mdremote.hiteqm, (mdremote.steps - mdremote.hiteqm), mdremote.extra_compute, mdremote.timestep, positive_diameter_in, negative_diameter_in);
                }
                else
                {
                  box.generate_lammps_datafile_chargedsurface(saltion_in, pz_in, nz_in, ion, smaller_ion_diameter, charge_meshpoint, counterions, valency_counterion, fraction_diameter, surface_area);
                  generateLammpsInputfileForChargedSurface(ein, mdremote.freq, mdremote.hiteqm, (mdremote.steps - mdremote.hiteqm), mdremote.extra_compute, mdremote.timestep, positive_diameter_in, negative_diameter_in);
                }

                cout << "Lammps Preprocessing ended." << endl;

            } else {

                cout << "Lammps Postprocessing started." << endl;

                int cnt_filename = 0;
                output_lammps(ion, cnt_filename, mdremote.freq);
                if (world.rank() == 0)
                    cout << "Number of samples used to get density profile:" << cnt_filename << endl;
                int lammps_samples = cnt_filename;
              //  make_bins(bin, box, bin_width);    // set up bins to be used for computing density profiles
                double lammps_density_profile_samples = 0;
                for (int cpmdstep = 0; cpmdstep < lammps_samples; cpmdstep++) {
                    vector <PARTICLE> ion;
                    vector<double> initial_density;
                    lammps_density_profile_samples++;
                    ReadParticlePositions(ion, cpmdstep, lammps_samples, positive_diameter_in, box, mdremote.freq);
                //    bin_ions(ion, box, initial_density, bin);
                    compute_density_profile(cpmdstep, lammps_density_profile_samples, mean_positiveion_density,
                                            mean_sq_positiveion_density,
                                            mean_negativeion_density, mean_sq_negativeion_density, ion, box, bin,
                                            mdremote, screen);

                }
                average_errorbars_density(lammps_density_profile_samples, mean_positiveion_density,
                                          mean_sq_positiveion_density,
                                          mean_negativeion_density,
                                          mean_sq_negativeion_density, ion, box, bin, simulationParams, screen);

                cout << "Lammps Postprocessing ended." << endl;

            }
        }
    }
    get_NetChargeDensity(simulationParams); // The net charge is zero if there is no charge on surfaces and the sizes of positive and negative ions are equal.
    if (charge_density != 0 && world.rank() == 0)
    {
      screen = true;
      cout << "Screen Factor Postprocessing started." << endl;
      int number_of_bins = int(box.lz / bin_width);
      bin_width = (box.lz / number_of_bins); // This give us the correct value of bin_width (see function.ccp)
      bin_width = bin_width * 0.01; // we choose smaller bin_width to calculate screen factor;
      int cnt_filename = 0;
      int samples = 0;
      bin.clear();
      make_bins(bin, box, bin_width);    // set up bins to be used for computing density profiles
      vector<double> mean_positiveion_density;            // average density profile
      vector<double> mean_negativeion_density;            // average density profile
      vector<double> mean_sq_positiveion_density;            // average of square of density
      vector<double> mean_sq_negativeion_density;            // average of square of density
      for (unsigned int b = 0; b < bin.size(); b++)
      {
          mean_positiveion_density.push_back(0.0);
          mean_negativeion_density.push_back(0.0);
          mean_sq_positiveion_density.push_back(0.0);
          mean_sq_negativeion_density.push_back(0.0);
      }
      double leftContact = -0.5 * box.lz + 0.5 * ion[0].diameter - 0.5 * bin[0].width;
      double rightContact = 0.5 * box.lz - 0.5 * ion[0].diameter - 0.5 * bin[0].width;
      bin[bin.size() - 1].lower = leftContact;
      bin[bin.size() - 2].lower = rightContact;
      bin[bin.size() - 1].higher = leftContact + bin[0].width;
      bin[bin.size() - 2].higher = rightContact + bin[0].width;
      bin[bin.size() - 1].midPoint = 0.5 * (bin[bin.size() - 1].lower + bin[bin.size() - 1].higher);
      bin[bin.size() - 2].midPoint = 0.5 * (bin[bin.size() - 2].lower + bin[bin.size() - 2].higher);

      output_lammps(ion, cnt_filename, mdremote.freq);  // if the simulation is done with lammps, this fuction rewrites the files in "temp" folder;
      if (world.rank() == 0)
          cout << "Number of samples used to get screen factor profile: " << cnt_filename << endl;
      samples = cnt_filename;
      double screen_density_profile_samples = 0;
      for (int cpmdstep = 0; cpmdstep < samples; cpmdstep++)
      {
        vector <PARTICLE> ion;
        vector<double> initial_density;
        screen_density_profile_samples++;
        ReadParticlePositions(ion, cpmdstep, samples, positive_diameter_in, box, mdremote.freq);
        bin_ions(ion, box, initial_density, bin);    // bin the ions to get initial density profile
        compute_density_profile(cpmdstep, screen_density_profile_samples, mean_positiveion_density,
                                mean_sq_positiveion_density,
                                mean_negativeion_density, mean_sq_negativeion_density, ion, box, bin,
                                mdremote, screen);
      }
      average_errorbars_density(screen_density_profile_samples, mean_positiveion_density,
                                mean_sq_positiveion_density,
                                mean_negativeion_density,
                                mean_sq_negativeion_density, ion, box, bin, simulationParams, screen);
      get_ScreeningFactor(charge_density, bin_width, simulationParams);
      cout << "Screen Factor Postprocessing ended." << endl;
    }



    return 0;
}

// End of main
