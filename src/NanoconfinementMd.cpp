// Last update: April 5, 2013
// This is main.
// This is MD simulation of an electrolyte confined within planar walls
// There is no dielectric contrast between the walls. This is uniform dielectric medium.
// Problem : Compute density profiles of ions trapped within planar walls
/* Useful studies :	
		     1. Role of valency of ions
		     2. Role of varying salt concentration
*/
#include "NanoconfinementMd.h"

int NanoconfinementMd::startSimulation(int argc, char *argv[], bool paraMap) {
    // Electrostatic system variables
    double bx, by, bz;        // lengths of the box
    double ein;            // permittivity of inside medium
    double eout;            // permittivity of outside medium
    int pz_in;            // positive valency of ions inside
    int nz_in;            // negative valency of ions inside
    double salt_conc_in;        // salt concentration outside	(enter in M)
    double saltion_diameter_in;    // inside salt ion diameter	(positive and negative ions assumed to have same diameter)
    double T;            // temperature at which the system of ions is

    // Simulation related variables
    double fraction_diameter;        // fraction that multiplies the diameter to generate the discretization width for the interface
    double Q;                // thermostat mass required to generate canonical ensemble
    unsigned int chain_length_real;    // Nose Hoover thermostat chain length for particles
    double bin_width;            // width of the bins used to compute density profiles
    CONTROL mdremote;            // remote control for md
    string config_file ="input_config.cfg";         //Configuration file path holder
    string simulationParams;
    char simulationParamsTemp[200];

    // Different parts of the system
    vector <PARTICLE> saltion_in;        // salt ions inside
    vector <PARTICLE> ion;            // all ions in the system
    INTERFACE box;            // interface, z planes hard walls; rest periodic boundaries

    // Analysis
    vector <DATABIN> bin;            // bins

    cout << "\nProgram starts\n";
    cout << "For help with the menu, type ./md_simulation_confined_ions -h" << endl;
    cout
            << "The default values are for a short run (~15 minutes on a single processor) that captures essential features; research-level simulations require larger number of steps (or compute time) which can be changed by changing -S. This leads to changing other parameters as well to avoid large-size data files."
            << endl;
    cout
            << "The default app simulates a total of 324 ions; for a faster run to test the app and understand what it produces select salt concentration (-c) to be 0.01 (so you can type and ENTER ./md_simulation_confined_ions -c 0.01). In this case, 32 salt ions corresponding to a density of 0.01 M are confined between two surfaces 3 nanometers apart. This should finish within 2 minutes on a single processor."
            << endl;
    cout
            << "OpenMP acceleration is built in the app which is essential for any research-intensive simulation; default is 1 processor."
            << endl;
    cout << "-----------------------------------------------------" << endl;

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
            ("epsilon_in,e", value<double>(&ein)->default_value(80.1),
             "dielectric const inside")        // must have ein = eout
            ("epsilon_out,E", value<double>(&eout)->default_value(80.1),
             "dielectric const outside")        // must have ein = eout
            ("fraction_diameter,g", value<double>(&fraction_diameter)->default_value(1),
             "for interface discretization width")    // enter a perfect square
            ("thermostat_mass,Q", value<double>(&Q)->default_value(1.0), "thermostat mass")
            ("chain_length_real,L", value<unsigned int>(&chain_length_real)->default_value(5),
             "chain length for real system: enter L+1 if you want L thermostats")
            ("bin_width,B", value<double>(&bin_width)->default_value(0.10), "bin width")
            ("md_timestep,T", value<double>(&mdremote.timestep)->default_value(0.0005), "time step used in md")
            ("md_eqm,P", value<int>(&mdremote.hiteqm)->default_value(100000), "production begin (md)")
            ("md_freq,F", value<int>(&mdremote.freq)->default_value(100), "sample frequency (md)")
            ("md_extra_compute,x", value<int>(&mdremote.extra_compute)->default_value(10000),
             "compute additional (md)")
            ("md_writedensity,w", value<int>(&mdremote.writedensity)->default_value(100000), "write density files")
            ("md_movie_freq,m", value<int>(&mdremote.moviefreq)->default_value(10000), "compute additional (md)")
            //("config,conf", value<string>(&config_file)->default_value("input_config.cfg"),
            // "name of a file of a configuration.")
            ;



    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    //-Z 3 -p 1 -n -1 -c 0.5 -d 0.714 -S 1000000
    options_description config("Configuration");
    config.add_options()
            ("confinement_length,Z", value<double>(&bz)->default_value(3),
             "box length in z direction in nanometers")        // enter in nanometers
            ("positive_valency,p", value<int>(&pz_in)->default_value(1), "positive valency inside")
            ("negative_valency,n", value<int>(&nz_in)->default_value(-1), "negative valency inside")
            ("salt_concentration,c", value<double>(&salt_conc_in)->default_value(0.50), "salt concentration inside")
            ("ion_diameter,d", value<double>(&saltion_diameter_in)->default_value(0.714),
             "salt ion diameter inside")        // enter in nanometers
            ("simulation_steps,S", value<int>(&mdremote.steps)->default_value(1000000), "steps used in md")
            ;

    options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    variables_map vm;
    store(parse_command_line(argc, argv, cmdline_options), vm);
    notify(vm);

    ifstream ifs(config_file.c_str());
    if (!ifs)
    {
        cout << "can not open config file: " << config_file << "\n";
        //return 0;
    }
    else
    {
        store(parse_config_file(ifs, config), vm);
        notify(vm);
    }

    if (vm.count("help")) {
        std::cout << cmdline_options << "\n";
        return 0;
    }

    //X and Y mapping
    if(paraMap){
    bx = sqrt(212/0.6022/salt_conc_in/bz);
    by=bx;
    mdremote.hiteqm=mdremote.steps*0.1;
    mdremote.writedensity=mdremote.steps*0.1;
    }

    //creating simulation parameters to create the output file names
    sprintf(simulationParamsTemp, "_%.2f_%d_%d_%.2f_%.3f_%d", bz,pz_in,nz_in,salt_conc_in,saltion_diameter_in,mdremote.steps);
    simulationParams=string(simulationParamsTemp);
    // Set up the system
    T = 1;        // set temperature
    box.ein = ein;
    box.eout = eout;
    box.set_up(salt_conc_in, 0, pz_in, 0, bx / unitlength, by / unitlength, bz / unitlength);
    box.put_saltions_inside(saltion_in, pz_in, nz_in, salt_conc_in, saltion_diameter_in, ion);
    make_bins(bin, box, bin_width);    // set up bins to be used for computing density profiles
    vector<double> initial_density;
    bin_ions(ion, box, initial_density, bin);    // bin the ions to get initial density profile

//   box.number_of_vertices = total_gridpoints;
    box.discretize(saltion_diameter_in / unitlength, fraction_diameter);

    // output to screen the parameters of the problem
    cout << "Dielectric constant of water is taken to be " << epsilon_water << endl;
#pragma omp parallel default(none)
    {
        if (omp_get_thread_num() == 0) {
            printf("Number of processors used (the app comes with OpenMP parallelization) %d\n",
                   omp_get_num_threads());
            printf("Make sure that number of grid points and ions is greater than %d\n", omp_get_num_threads());
        }
    }
    //cout << "Number of processors used (the app comes with OpenMP parallelization)  " << THREADSIZE << endl;
    //cout << "Make sure that number of grid points and ions is greater than  " << THREADSIZE << endl;
    cout << "Unit of length is " << unitlength << " nanometers"
         << endl; // half Bjerrum length; close to Na ion radius
    cout << "Unit of mass is " << unitmass << " grams" << endl; // mass of sodium atom in CGS; in grams
    cout << "Unit of energy is " << unitenergy << " ergs (CGS)" << endl; // kB room_T
    cout << "Unit of time is " << unittime << " s"
         << endl; // reduced units (LJ); you can derive from above 3 units.
    cout << "Reduced units: scalefactor entering in Coloumb interaction is " << scalefactor << endl;
    cout << "Box dimensions (in reduced units) x | y | z " << setw(15) << box.lx << setw(15) << box.ly << setw(15)
         << box.lz << endl;
    cout << "Box dimensions (in nanometers) x | y | z " << setw(15) << bx << setw(15) << by << setw(15) << bz
         << endl;
    cout << "Permittivity inside " << box.ein << endl;
    cout << "Permittivity outside " << box.eout << endl;
    cout << "Contrast strength " << 2 * (box.eout - box.ein) / (box.eout + box.ein) << endl;
    cout << "Positive ion valency inside " << pz_in << endl;
    cout << "Negative ion valency inside " << nz_in << endl;
    cout << "Salt ion diameter inside " << saltion_diameter_in / unitlength << endl;
    cout << "Salt concentration inside " << salt_conc_in << " M" << endl;
    cout << "Debye length inside " << box.inv_kappa_in << endl;
    cout << "Mean separation inside " << box.mean_sep_in << endl;
    cout << "Number of salt ions inside " << saltion_in.size() << endl;
    cout << "Temperature " << T << endl;
    cout << "Binning width (uniform) " << bin[0].width << endl;
    cout << "Number of bins " << bin.size() << endl;
    cout << "Number of points discretizing the left and right z planar walls/interfaces/surfaces "
         << box.leftplane.size() << "  " << box.rightplane.size() << endl;

    // write to files

    // initial configuration
	string initial_configurationPath= rootDirectory+"outfiles/initialconfig.dat";
    ofstream initial_configuration(initial_configurationPath.c_str());
    for (unsigned int i = 0; i < ion.size(); i++)
        initial_configuration << "ion" << setw(5) << ion[i].id << setw(15) << "charge" << setw(5) << ion[i].q
                              << setw(15) << "position" << setw(15) << ion[i].posvec << endl;
    initial_configuration.close();

    // initial density
	string density_profilePath= rootDirectory+"outfiles/initial_density_profile.dat";
    ofstream density_profile(density_profilePath.c_str(), ios::out);
    for (unsigned int b = 0; b < initial_density.size(); b++)
        density_profile << bin[b].lower << setw(15) << initial_density.at(b) << endl;
    density_profile.close();

    // check point
    double totalions = 0;
    for (unsigned int b = 0; b < initial_density.size(); b++)
        totalions += initial_density.at(b) * bin[b].volume;
    cout << "total ions " << totalions << ", should be " << ion.size() << endl;
    int totalpions = 0, totalnions = 0;
    for (unsigned int i = 0; i < ion.size(); i++) {
        if (ion[i].valency > 0)
            totalpions += 1;
        else if (ion[i].valency < 0)
            totalnions += 1;
    }
    cout << "Total positive ions " << totalpions << endl;
    cout << "Total negative ions " << totalnions << endl;

    // some calculations before simulation begins
    cout << "Total charge inside the sphere " << box.total_charge_inside(ion) << endl;

    // prepare for md : make real baths
    vector <THERMOSTAT> real_bath;
    if (chain_length_real == 1)
        real_bath.push_back(THERMOSTAT(0, T, 3 * ion.size(), 0.0, 0, 0));
    else {
        real_bath.push_back(THERMOSTAT(Q, T, 3 * ion.size(), 0, 0, 0));
        while (real_bath.size() != chain_length_real - 1)
            real_bath.push_back(THERMOSTAT(Q / (3 * ion.size()), T, 1, 0, 0, 0));
        real_bath.push_back(THERMOSTAT(0, T, 3 * ion.size(), 0.0, 0, 0));
        // final bath is dummy bath (dummy bath always has zero mass)
    }

    // Car-Parrinello Molecular Dynamics
    md(ion, box, real_bath, bin, mdremote,simulationParams);

    // Post simulation analysis (useful for short runs, but performed otherwise too)
    cout << "MD trust factor R (should be < 0.05) is " << compute_MD_trust_factor_R(mdremote.hiteqm) << endl;
    // Perform the following calculation when testing for how frequently you should sample data to ensure the samples are decorrelated
    //auto_correlation_function();

    cout << "Program ends \n\n";

    return 0;
}

// End of main