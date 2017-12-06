// This file contains the routine that computes the force on the particle i for all i

#include "forces.h"

// Forces on particles
void for_md_calculate_force(vector <PARTICLE> &ion, INTERFACE &box, char flag, vector <vector<VECTOR3D> > &forceMatrix, vector<VECTOR3D> &electrostatic_forces, vector<VECTOR3D> &lj_ion_ion, vector<VECTOR3D> &lj_ion_leftdummy, vector<VECTOR3D> &lj_ion_left_wall, vector<VECTOR3D> &lj_ion_rightdummy, vector<VECTOR3D> &lj_ion_right_wall) {

    unsigned int i, j, k, iloop, j1;
    int factor;
    VECTOR3D h1, temp_vec;
    double hcsh, E_z, r1, r2, dz;
    long double r, r3;
    double elj, r6, r12, d, d2, d6, d12;

    // force on the particles (electrostatic)
    // parallel calculation of forces (uniform case)
#pragma omp parallel default(shared) private(iloop, j1, h1, dz, factor, r1, r2, E_z, hcsh, temp_vec, r, r3)
    {
#pragma omp for schedule(dynamic) nowait
        for (iloop = 0; iloop < ion.size(); iloop++) {
            for (j1 = iloop+1; j1 < ion.size(); j1++) {
                h1 = VECTOR3D(0, 0, 0);
                dz = ion[iloop].posvec.z - ion[j1].posvec.z;
                if (dz >= 0) factor = 1;
                else factor = -1;
                r1 = sqrt(0.5 + (dz / box.lx) * (dz / box.lx));
                r2 = sqrt(0.25 + (dz / box.lx) * (dz / box.lx));
                E_z = 4 * atan(4 * fabs(dz) * r1 / box.lx);
                hcsh = (4 / box.lx) * (1 / (r1 * (0.5 + r1)) - 1 / (r2 * r2)) * dz + factor * E_z +
                       16 * fabs(dz) * (box.lx / (box.lx * box.lx + 4 * dz * dz * r1 * r1)) *
                       (fabs(dz) * dz / (box.lx * box.lx * r1) + factor * r1);

                h1.z = h1.z + 2 * ion[iloop].q * (ion[j1].q / (box.lx * box.lx)) * 0.5 *
                              (1 / ion[iloop].epsilon + 1 / ion[j1].epsilon) * hcsh;

                temp_vec = ion[iloop].posvec - ion[j1].posvec;
                if (temp_vec.x > box.lx / 2) temp_vec.x -= box.lx;
                if (temp_vec.x < -box.lx / 2) temp_vec.x += box.lx;
                if (temp_vec.y > box.ly / 2) temp_vec.y -= box.ly;
                if (temp_vec.y < -box.ly / 2) temp_vec.y += box.ly;
                r = temp_vec.GetMagnitude();
                r3 = r * r * r;
                h1 = h1 + ((temp_vec ^ ((-1.0) / r3)) ^
                           ((-0.5) * ion[iloop].q * ion[j1].q * (1 / ion[iloop].epsilon + 1 / ion[j1].epsilon)));
                // force is q1 * q2 / r^2, no half factor. half factor cancels out in the above expression.
                forceMatrix[iloop][j1] = h1;
            }

        }
    }

#pragma omp parallel for schedule(dynamic) private(iloop, j1, h1)
    for (iloop = 0; iloop < ion.size(); iloop++) {
        h1 = VECTOR3D(0, 0, 0);
        for (j1 = 0; j1 < ion.size(); j1++) {
            if (j1 == iloop) continue;
            if (iloop<j1)
                h1 = h1 + forceMatrix[iloop][j1];
            else
                h1 = h1 + (forceMatrix[j1][iloop]^ (-1));
        }
        electrostatic_forces[iloop]=h1^(scalefactor);
    }

    // excluded volume interactions given by purely repulsive LJ
    // ion-ion
    VECTOR3D r_vec, flj, fljcc;
    PARTICLE wall_dummy, dummy;
#pragma omp parallel default(shared) private(i, j, r_vec, r2, d, d2, elj, r6, r12, d6, d12, fljcc)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = 0; i < ion.size(); i++) {
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
            lj_ion_ion[i] = fljcc;
        }
    }
    // ion-box

    // interaction with the left plane hard wall

    // make a dummy particle with the same diameter as the ion and touching left of the left wall s. t. it is closest to the ion
    for (unsigned int i = 0; i < ion.size(); i++) {
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
        lj_ion_leftdummy[i] = flj;
    }

    // ion interacting with discretized left wall
#pragma omp parallel default(shared) private(i, k, wall_dummy, r_vec, r2, d, d2, elj, r6, r12, d6, d12, flj)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = 0; i < ion.size(); i++) {
            flj = VECTOR3D(0, 0, 0);
            if (ion[i].posvec.z < -0.5 * box.lz +
                                  ion[i].diameter)  // avoiding calculating interactions between left wall and ions in bulk. -Yufei - Vikram
            {
                for (k = 0; k < box.leftplane.size(); k++) {
                    wall_dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.eout,
                                          VECTOR3D(box.leftplane[k].posvec.x, box.leftplane[k].posvec.y,
                                                   box.leftplane[k].posvec.z - 0.5 * ion[i].diameter), box.lx, box.ly,
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
            lj_ion_left_wall[i] = flj;
        }
    }

    //interaction with the right plane hard wall

    //make a dummy particle with the same diameter as the ion and touching right of the right wall s. t. it is closest to the ion
    for (unsigned int i = 0; i < ion.size(); i++) {
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
        lj_ion_rightdummy[i] = flj;
    }

    // ion interacting with discretized right wall
#pragma omp parallel default(shared) private(i, k, wall_dummy, r_vec, r2, d, d2, elj, r6, r12, d6, d12, flj)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = 0; i < ion.size(); i++) {
            flj = VECTOR3D(0, 0, 0);
            if (ion[i].posvec.z > 0.5 * box.lz -
                                  ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk. -Yufei -Vikram
            {
                for (k = 0; k < box.rightplane.size(); k++) {
                    wall_dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.eout,
                                          VECTOR3D(box.rightplane[k].posvec.x, box.rightplane[k].posvec.y,
                                                   box.rightplane[k].posvec.z + 0.5 * ion[i].diameter),
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
            lj_ion_right_wall[i] = flj;
        }
    }

    // total force on the particle = the electrostatic force + the Lennard-Jones force
    //h1 = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < ion.size(); i++)
        ion[i].forvec =
                electrostatic_forces[i]  + lj_ion_ion[i] + lj_ion_leftdummy[i] + lj_ion_left_wall[i] + lj_ion_rightdummy[i] +
                lj_ion_right_wall[i];
        //electrostatic_forces[i]=h1;


    return;
}