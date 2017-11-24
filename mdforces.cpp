// This file contains the routine that computes the force on the particle i for all i

#include "forces.h"

// Forces on particles
void for_md_calculate_force(vector <PARTICLE> &ion, INTERFACE &box, char flag, unsigned int lowerBound,
                            unsigned int upperBound, vector <VECTOR3D> &partialForceVector,
                            vector <VECTOR3D> &lj_ion_ion, vector <VECTOR3D> &lj_ion_leftdummy,
                            vector <VECTOR3D> &lj_ion_left_wall, vector <VECTOR3D> &lj_ion_rightdummy,
                            vector <VECTOR3D> &lj_ion_right_wall, vector <VECTOR3D> &sendForceVector) {

    mpi::environment env;
    mpi::communicator world;

    unsigned int i, j;
    int factor;
    VECTOR3D h1, temp_vec;
    double hcsh, E_z, r1, r2, dz;
    long double r, r3;
    double elj, r6, r12, d, d2, d6, d12;
    // force on the particles (electrostatic)
    // parallel calculation of forces (uniform case)

#pragma omp parallel default(shared) private(i, j, h1, dz, factor, r1, r2, E_z, hcsh, temp_vec, r, r3)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = lowerBound; i <= upperBound; i++) {
            h1 = VECTOR3D(0, 0, 0);
            for (j = 0; j < ion.size(); j++) {
                if (j == i) continue;
                dz = ion[i].posvec.z - ion[j].posvec.z;
                if (dz >= 0) factor = 1;
                else factor = -1;
                r1 = sqrt(0.5 + (dz / box.lx) * (dz / box.lx));
                r2 = sqrt(0.25 + (dz / box.lx) * (dz / box.lx));
                E_z = 4 * atan(4 * fabs(dz) * r1 / box.lx);
                hcsh = (4 / box.lx) * (1 / (r1 * (0.5 + r1)) - 1 / (r2 * r2)) * dz + factor * E_z +
                       16 * fabs(dz) * (box.lx / (box.lx * box.lx + 4 * dz * dz * r1 * r1)) *
                       (fabs(dz) * dz / (box.lx * box.lx * r1) + factor * r1);

                h1.z = h1.z +
                       2 * ion[i].q * (ion[j].q / (box.lx * box.lx)) * 0.5 * (1 / ion[i].epsilon + 1 / ion[j].epsilon) *
                       hcsh;

                temp_vec = ion[i].posvec - ion[j].posvec;
                if (temp_vec.x > box.lx / 2) temp_vec.x -= box.lx;
                if (temp_vec.x < -box.lx / 2) temp_vec.x += box.lx;
                if (temp_vec.y > box.ly / 2) temp_vec.y -= box.ly;
                if (temp_vec.y < -box.ly / 2) temp_vec.y += box.ly;
                r = temp_vec.GetMagnitude();
                r3 = r * r * r;
                h1 =h1+ ((temp_vec ^ ((-1.0) / r3)) ^
                      ((-0.5) * ion[i].q * ion[j].q * (1 / ion[i].epsilon + 1 / ion[j].epsilon)));

                // force is q1 * q2 / r^2, no half factor. half factor cancels out in the above expression.
            }
            sendForceVector[i - lowerBound] = ((h1) ^ (scalefactor));


        }
    }

    // excluded volume interactions given by purely repulsive LJ
    // ion-ion
    //vector<VECTOR3D> lj_ion_ion (ion.size(), VECTOR3D(0,0,0));
    VECTOR3D r_vec, flj, fljcc;
    PARTICLE wall_dummy, dummy;
#pragma omp parallel default(shared) private(i, j, r_vec, r2, d, d2, elj, r6, r12, d6, d12, fljcc)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = lowerBound; i <= upperBound; i++) {
            fljcc = VECTOR3D(0, 0, 0);
            for (j = 0; j < ion.size(); j++) {
                if (j == i) continue;
                r_vec = ion[i].posvec - ion[j].posvec;
                if (r_vec.x > box.lx / 2) r_vec.x -= box.lx;
                if (r_vec.x < -box.lx / 2) r_vec.x += box.lx;
                if (r_vec.y > box.ly / 2) r_vec.y -= box.ly;
                if (r_vec.y < -box.ly / 2) r_vec.y += box.ly;
                r2 = r_vec.GetMagnitudeSquared();
                d = 0.5 * (ion[i].diameter + ion[j].diameter);
                d2 = d * d;
                elj = 1.0;
                if (r2 < dcut2 * d2) {
                    r6 = r2 * r2 * r2;
                    r12 = r6 * r6;
                    d6 = d2 * d2 * d2;
                    d12 = d6 * d6;
                    fljcc = fljcc + (r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
                } else
                    fljcc = fljcc + VECTOR3D(0, 0, 0);
            }
            lj_ion_ion[i - lowerBound] = fljcc;
        }
    }
    // ion-box
    // interaction with the left plane hard wall

    // make a dummy particle with the same diameter as the ion and touching left of the left wall s. t. it is closest to the ion
    //vector<VECTOR3D> lj_ion_leftdummy (ion.size(), VECTOR3D(0,0,0));

//#pragma omp parallel for schedule(dynamic) private(i, r_vec, r2, d, d2, elj, r6, r12, d6, d12, flj)
    for (i = lowerBound; i <= upperBound; i++) {
        flj = VECTOR3D(0, 0, 0);
        if (ion[i].posvec.z < -0.5 * box.lz +
                              ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk. replacing 1 by diameter. -Yufei -Vikram -Vikram
        {
            dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.eout,
                             VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, -0.5 * box.lz - 0.5 * ion[i].diameter), box.lx,
                             box.ly, box.lz);
            r_vec = ion[i].posvec - dummy.posvec;
            r2 = r_vec.GetMagnitudeSquared();
            d = 0.5 * (ion[i].diameter + (dummy.diameter));
            d2 = d * d;
            elj = 1.0;
            if (r2 < dcut2 * d2) {
                r6 = r2 * r2 * r2;
                r12 = r6 * r6;
                d6 = d2 * d2 * d2;
                d12 = d6 * d6;
                flj = r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
            }
        }
        lj_ion_leftdummy[i - lowerBound] = flj;

    }

    // ion interacting with discretized left wall
    //vector<VECTOR3D> lj_ion_left_wall (ion.size(), VECTOR3D(0,0,0));
#pragma omp parallel default(shared) private(i, j, wall_dummy, r_vec, r2, d, d2, elj, r6, r12, d6, d12, flj)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = lowerBound; i <= upperBound; i++) {
            flj = VECTOR3D(0, 0, 0);
            if (ion[i].posvec.z < -0.5 * box.lz +
                                  ion[i].diameter)  // avoiding calculating interactions between left wall and ions in bulk. -Yufei - Vikram
            {
                for (j = 0; j < box.leftplane.size(); j++) {
                    wall_dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.eout,
                                          VECTOR3D(box.leftplane[j].posvec.x, box.leftplane[j].posvec.y,
                                                   box.leftplane[j].posvec.z - 0.5 * ion[i].diameter), box.lx, box.ly,
                                          box.lz);
                    r_vec = ion[i].posvec - wall_dummy.posvec;
                    if (r_vec.x > box.lx / 2) r_vec.x -= box.lx;
                    if (r_vec.x < -box.lx / 2) r_vec.x += box.lx;
                    if (r_vec.y > box.ly / 2) r_vec.y -= box.ly;
                    if (r_vec.y < -box.ly / 2) r_vec.y += box.ly;
                    r2 = r_vec.GetMagnitudeSquared();
                    d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
                    d2 = d * d;
                    elj = 1.0;
                    if (r2 < dcut2 * d2) {
                        r6 = r2 * r2 * r2;
                        r12 = r6 * r6;
                        d6 = d2 * d2 * d2;
                        d12 = d6 * d6;
                        flj += r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
                    }
                }
            }
            lj_ion_left_wall[i - lowerBound] = flj;
        }
    }
    //interaction with the right plane hard wall

    //make a dummy particle with the same diameter as the ion and touching right of the right wall s. t. it is closest to the ion
    //vector<VECTOR3D> lj_ion_rightdummy (ion.size(), VECTOR3D(0,0,0));
    //#pragma omp parallel for schedule(dynamic) private(i, r_vec, r2, d, d2, elj, r6, r12, d6, d12, flj)
    for (i = lowerBound; i <= upperBound; i++) {
        flj = VECTOR3D(0, 0, 0);
        if (ion[i].posvec.z > 0.5 * box.lz -
                              ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
        {
            dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.eout,
                             VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, 0.5 * box.lz + 0.5 * ion[i].diameter), box.lx,
                             box.ly, box.lz);
            r_vec = ion[i].posvec - dummy.posvec;
            r2 = r_vec.GetMagnitudeSquared();
            d = 0.5 * (ion[i].diameter + dummy.diameter);
            d2 = d * d;
            elj = 1.0;
            if (r2 < dcut2 * d2) {
                r6 = r2 * r2 * r2;
                r12 = r6 * r6;
                d6 = d2 * d2 * d2;
                d12 = d6 * d6;
                flj = r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
            }
        }
        lj_ion_rightdummy[i - lowerBound] = flj;

    }

    // ion interacting with discretized right wall
//vector<VECTOR3D> lj_ion_right_wall (ion.size(), VECTOR3D(0,0,0));
#pragma omp parallel default(shared) private(i, j, wall_dummy, r_vec, r2, d, d2, elj, r6, r12, d6, d12, flj)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = lowerBound; i <= upperBound; i++) {
            flj = VECTOR3D(0, 0, 0);
            if (ion[i].posvec.z > 0.5 * box.lz -
                                  ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
            {
                for (j = 0; j < box.rightplane.size(); j++) {
                    wall_dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.eout,
                                          VECTOR3D(box.rightplane[j].posvec.x, box.rightplane[j].posvec.y,
                                                   box.rightplane[j].posvec.z + 0.5 * ion[i].diameter),
                                          box.lx, box.ly, box.lz);
                    r_vec = ion[i].posvec - wall_dummy.posvec;
                    if (r_vec.x > box.lx / 2) r_vec.x -= box.lx;
                    if (r_vec.x < -box.lx / 2) r_vec.x += box.lx;
                    if (r_vec.y > box.ly / 2) r_vec.y -= box.ly;
                    if (r_vec.y < -box.ly / 2) r_vec.y += box.ly;
                    r2 = r_vec.GetMagnitudeSquared();
                    d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
                    d2 = d * d;
                    elj = 1.0;
                    if (r2 < dcut2 * d2) {
                        r6 = r2 * r2 * r2;
                        r12 = r6 * r6;
                        d6 = d2 * d2 * d2;
                        d12 = d6 * d6;
                        flj += r_vec ^ (48 * elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2));
                    }
                }
            }
            lj_ion_right_wall[i - lowerBound] = flj;
        }
    }


    if (world.size() > 1) {

        // total force on the particle = the electrostatic force + the Lennard-Jones force slave processes
        for (i = 0; i < sendForceVector.size(); i++) {
            sendForceVector[i] = sendForceVector[i] + lj_ion_ion[i] + lj_ion_leftdummy[i] +
                                 lj_ion_left_wall[i] + lj_ion_rightdummy[i] + lj_ion_right_wall[i];
        }

        //broadcasting using all gather = gather + broadcast
        //template<typename T> void all_gather(const communicator & comm, const T * in_values, int n, std::vector< T > & out_values);
         all_gather(world, &sendForceVector[0],sendForceVector.size(),partialForceVector);

        //assigning all gathered forces to each node's force vector
        //#pragma omp parallel for schedule(dynamic) private(i)
        for (i = 0; i < ion.size(); i++)
            ion[i].forvec = partialForceVector[i];

    } else {

        // total force on the particle = the electrostatic force + the Lennard-Jones force in main processes
        for (i = lowerBound; i <= upperBound; i++)
            ion[i].forvec =
                    sendForceVector[i - lowerBound] + lj_ion_ion[i - lowerBound] + lj_ion_leftdummy[i - lowerBound] +
                    lj_ion_left_wall[i - lowerBound] + lj_ion_rightdummy[i - lowerBound] +
                    lj_ion_right_wall[i - lowerBound];

    }

    return;
}