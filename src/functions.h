// This is a header file containing functions

#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include "utility.h"
#include "interface.h"
#include "particle.h"
#include "vertex.h"
#include "control.h"
#include "databin.h"
#include "thermostat.h"
#include "energies.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// general functions
// -----------------

// display progress bar (code from the internet)
void ProgressBar(double);

// overloaded << to print 3d vectors
ostream &operator<<(ostream &, VECTOR3D);

// make bins
void make_bins(vector<DATABIN> &, INTERFACE &, double);

// bin ions
void bin_ions(vector<PARTICLE> &, INTERFACE &, vector<double> &, vector<DATABIN> &);

// initialize particle velocities
void initialize_particle_velocities(vector<PARTICLE> &, vector<THERMOSTAT> &);

// compute and write useful data in cpmd
void compute_n_write_useful_data(int, vector<PARTICLE> &, vector<THERMOSTAT> &, INTERFACE &, unsigned int, unsigned int,
                                 vector<double> &, vector<double> &,
                                 vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector <double> &, vector <double> &, double, int);

// make movie
void make_movie(int num, vector<PARTICLE> &ion, INTERFACE &box);

// compute density profile
void compute_density_profile(int cpmdstep, double density_profile_samples,
                             vector<double> &meanPositiveionDensity,
                             vector<double> &mean_sq_positiveion_density,
                             vector<double> &meanNegativeionDensity,
                             vector<double> &mean_sq_negativeion_density,
                             vector<PARTICLE> &ion, INTERFACE &box,
                             vector<DATABIN> &bin, CONTROL &cpmdremote, bool);

void average_errorbars_density(double density_profile_samples, vector<double>& mean_positiveion_density,
                                   vector<double>& mean_sq_positiveion_density,
                                   vector<double>& mean_negativeion_density,
                                   vector<double>& mean_sq_negativeion_density,
                                   vector<PARTICLE>& ion, INTERFACE &box,
                                   vector<DATABIN>& bin,  string simulationParams, bool);

void ReadParticlePositions(vector<PARTICLE>&, int, int, double, INTERFACE &, double);
void output_lammps(vector<PARTICLE>&, int&, double);
// post analysis : compute R
double compute_MD_trust_factor_R(int);

// post analysis : auto correlation function
void auto_correlation_function();

// generate LAMMPS input script
void generateLammpsInputfileForChargedSurface(double , int , int , int , int, double, double, double);
void generateLammpsInputfileForUnchargedSurface(double , int , int , int , int, double, double, double);
void get_NetChargeDensity(string);
void get_ScreeningFactor(double, double, string);

// functions useful in computing forces and energies
// -------------------------------------------------

// Green's function G
inline long double G(vector<VERTEX> &s, unsigned int k, unsigned int l) {
    if (l != k) return 1.0 / (s[k].posvec - s[l].posvec).GetMagnitude();
    if (l == k) return 2 * sqrt(pi) * sqrt(s[k].a) / s[k].a;
    return 0.0;
}

// computes gradient of green's fn
inline VECTOR3D Grad(VECTOR3D &vec1, VECTOR3D &vec2) {
    long double r = (vec1 - vec2).GetMagnitude();
    long double r3 = r * r * r;
    return (vec1 - vec2) ^ ((-1.0) / r3);
}

// Gradient of Green's function denoted by H
inline long double H(vector<VERTEX> &s, unsigned int k, unsigned int l, double radius) {
    if (l != k) return s[l].normalvec * Grad(s[l].posvec, s[k].posvec);
    if (l == k) return -0.5 * sqrt(pi) * sqrt(s[k].a) / (radius * s[k].a);
    return 0.0;
}

// computes gradient of normal dot gradient of 1/r
inline VECTOR3D GradndotGrad(VECTOR3D &vec1, VECTOR3D &vec2, VECTOR3D &normal) {
    long double r = (vec1 - vec2).GetMagnitude();
    long double r3 = r * r * r;
    long double r5 = r3 * r * r;
    return ((normal ^ (1.0 / r3)) - ((vec1 - vec2) ^ (3 * (normal * (vec1 - vec2)) / r5)));
}

// functions useful in implementing constraint
// -------------------------------------------

// constraint equation
inline long double constraint(vector<VERTEX> &s, vector<PARTICLE> &ion, INTERFACE &dsphere) {
    return dsphere.total_induced_charge(s) - dsphere.total_charge_inside(ion) * (1 / dsphere.eout - 1 / dsphere.ein);
}

// SHAKE to ensure constraint is true
inline void SHAKE(vector<VERTEX> &s, vector<PARTICLE> &ion, INTERFACE &dsphere,
                  CONTROL &simremote)    // remote of the considered simulation
{
    long double sigma = constraint(s, ion, dsphere);
    for (unsigned int k = 0; k < s.size(); k++)
        s[k].vw = s[k].vw - (1.0 / simremote.timestep) * sigma / (s[k].a * int(s.size()));
    for (unsigned int k = 0; k < s.size(); k++)
        s[k].w = s[k].w - sigma / (s[k].a * int(s.size()));
    return;
}

// dot constraint equation
inline long double dotconstraint(vector<VERTEX> &s) {
    long double sigmadot = 0;
    for (unsigned int k = 0; k < s.size(); k++)
        sigmadot += s[k].vw * s[k].a;
    return sigmadot;
}

// RATTLE to ensure time derivative of the constraint is true
inline void RATTLE(vector<VERTEX> &s) {
    long double sigmadot = dotconstraint(s);
    for (unsigned int k = 0; k < s.size(); k++)
        s[k].vw = s[k].vw - sigmadot / (s[k].a * int(s.size()));
    return;
}

// functions useful in Nose-Hoover chain implementation
// -------------------------------------------

// update bath xi value
inline void update_chain_xi(unsigned int j, vector<THERMOSTAT> &bath, double dt, long double ke) {
    if (bath[j].Q == 0)
        return;
    if (j != 0)
        bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j + 1].xi) + 0.5 * dt * (1.0 / bath[j].Q) *
                                                                    (bath[j - 1].Q * bath[j - 1].xi * bath[j - 1].xi -
                                                                     bath[j].dof * kB * bath[j].T) *
                                                                    exp(-0.25 * dt * bath[j + 1].xi);
//     bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j+1].xi) + 0.5 * dt * (1.0 / bath[j].Q) * (bath[j-1].Q * bath[j-1].xi * bath[j-1].xi - bath[j].dof * kB * bath[j].T) * (exp(-0.5 * dt * bath[j+1].xi) - 1) / (-0.5 * dt * bath[j+1].xi);
    else
        bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j + 1].xi) +
                     0.5 * dt * (1.0 / bath[j].Q) * (2 * ke - bath[j].dof * kB * bath[j].T) *
                     exp(-0.25 * dt * bath[j + 1].xi);
//     bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j+1].xi) + 0.5 * dt * (1.0 / bath[j].Q) * (2*ke - bath[j].dof * kB * bath[j].T) * (exp(-0.5 * dt * bath[j+1].xi) - 1) / (-0.5 * dt * bath[j+1].xi);
    return;
}

// functions that write data files
// -------------------------------

// bin ions to get density profile
inline void bin_ions(vector<PARTICLE> &ion, INTERFACE &box, vector<double> &density, vector<DATABIN> &bin) {
    double r;
    int bin_number;

    for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
        bin[bin_num].n = 0;
    for (unsigned int i = 0; i < ion.size(); i++) {
        /*0th bin is for left wall; bin 0 starts at -0.5*box.lz respect to ion origin */
        r = 0.5 * box.lz +
            ion[i].posvec.z;    // r should be positive, counting from half lz and that is the adjustment here; bin 0 corresponds to left wall
        bin_number = int(r / bin[0].width);
        bin[bin_number].n = bin[bin_number].n + 1;
    }
    /*This is to get contact point densities*/
    for (unsigned int i = 0; i < ion.size(); i++) {
        /*0th bin is for left wall; bin 0 starts at -0.5*box.lz respect to ion origin */
        r = 0.5 * box.lz + ion[i].posvec.z;
        if (bin[bin.size() - 2].lower <= ion[i].posvec.z && ion[i].posvec.z < bin[bin.size() - 2].higher)
            bin[bin.size() - 2].n = bin[bin.size() - 2].n + 1;
        else if (bin[bin.size() - 1].lower <= ion[i].posvec.z && ion[i].posvec.z < bin[bin.size() - 1].higher)
            bin[bin.size() - 1].n = bin[bin.size() - 1].n + 1;

    }

    for (unsigned int bin_num = 0; bin_num < bin.size(); bin_num++)
        density.push_back(
                bin[bin_num].n / bin[bin_num].volume);            // push_back is the culprit, array goes out of bound
    // volume now measured in inverse Molars, such that density is in M
    return;
}

#endif
