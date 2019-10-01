// This file contains the routine to compute the total potential energy of the system
// Electrostatic and Excluded volume (LJ) contributions

#include "energies.h"

// Potential energy
long double energy_functional(vector <PARTICLE> &ion, INTERFACE &box, unsigned int lowerBound,
                              unsigned int upperBound, vector <double> &ion_energy,
                              vector <double> &lj_ion_ion, vector <double> &lj_ion_leftdummy,
                              vector <double> &lj_ion_leftwall, vector <double> &lj_ion_rightdummy,
                              vector <double> &lj_ion_rightwall, vector <double> &coulomb_rightwall, vector <double> &coulomb_leftwall, double charge_meshpoint, int valency_counterion) {

	mpi::environment env;
	mpi::communicator world;
    unsigned int i, j, k, m, l, b;
    double fqq, fqq_csh, fqq_rightwall, fqq_leftwall, fqq_csh_leftwall,fqq_csh_rightwall;
    double fcsh_z, r1, r2, dz, fcsh_inf, dz_rightwall, r1_rightwall, r2_rightwall,dz_leftwall, r1_leftwall, fcsh_inf_rightwall, fcsh_z_rightwall, r2_leftwall, fcsh_z_leftwall, fcsh_inf_leftwall;
    VECTOR3D temp_vec, temp_vec_rightwall, temp_vec_leftwall;
    double potential,totalPotential;

    // push_back is not compatible with pragma, so we have initialize vector in this way.

    // charged sheets method is used to compute Coulomb interactions; an Ewald version should be designed to compare and ensure that long-range effects are taken into account in either methods

#pragma omp parallel default(shared) private(i, j, fqq, fqq_csh, dz, r1, r2, fcsh_inf, fcsh_z, temp_vec)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = lowerBound; i <= upperBound; i++) {
            fqq = 0;
            fqq_csh = 0;
            for (j = 0; j < ion.size(); j++) {
                dz = ion[i].posvec.z - ion[j].posvec.z;
                r1 = sqrt(0.5 + (dz / box.lx) * (dz / box.lx));
                r2 = sqrt(0.25 + (dz / box.lx) * (dz / box.lx));
                fcsh_z = 4 * box.lx * log((0.5 + r1) / r2) - fabs(dz) * (2 * pi - 4 * atan(4 * fabs(dz) * r1 / box.lx));
                fcsh_inf = -2 * pi * fabs(dz);
                fqq_csh += ion[i].q * (ion[j].q / (box.lx * box.lx)) * 0.5 * (1 / ion[i].epsilon + 1 / ion[j].epsilon) *
                           (fcsh_inf - fcsh_z);

                if (i == j) continue;

                temp_vec = ion[i].posvec - ion[j].posvec;

                if (temp_vec.x > box.lx / 2) temp_vec.x -= box.lx;
                if (temp_vec.x < -box.lx / 2) temp_vec.x += box.lx;
                if (temp_vec.y > box.ly / 2) temp_vec.y -= box.ly;
                if (temp_vec.y < -box.ly / 2) temp_vec.y += box.ly;

                fqq += 0.5 * ion[i].q * ion[j].q * 0.5 * (1.0 / ion[i].epsilon + 1.0 / ion[j].epsilon) /
                       ((temp_vec).GetMagnitude());
            }
            ion_energy[i-lowerBound] = fqq + fqq_csh;
        }
    }

    // Excluded volume interaction energy given by purely repulsive LJ

    // ion-ion
    VECTOR3D r_vec;
    PARTICLE wall_dummy;
    double uljcc, elj, r, r6, d, d2, d6;
#pragma omp parallel default(shared) private(i, j, uljcc, r, r2, r6, d, d2, d6, elj, r_vec)
    {
#pragma omp for schedule(dynamic) nowait
        for (i = lowerBound; i <= upperBound; i++) {
            uljcc = 0.0;
            for (j = 0; j < ion.size(); j++) {
                if (j == i) continue;
                r_vec = ion[i].posvec - ion[j].posvec;

                if (r_vec.x > box.lx / 2) r_vec.x -= box.lx;
                if (r_vec.x < -box.lx / 2) r_vec.x += box.lx;
                if (r_vec.y > box.ly / 2) r_vec.y -= box.ly;
                if (r_vec.y < -box.ly / 2) r_vec.y += box.ly;

                r = r_vec.GetMagnitude();
                d = 0.5 * (ion[i].diameter + ion[j].diameter);
                elj = 1.0;
                if (r < dcut * d) {
                    r2 = r * r;
                    r6 = r2 * r2 * r2;
                    d2 = d * d;
                    d6 = d2 * d2 * d2;
                    uljcc = uljcc + 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
                } else
                    uljcc = uljcc + 0.0;
            }
            lj_ion_ion[i-lowerBound] = uljcc;
        }
    }
    // ion-box

    // left wall

    // ion interacting with left wall directly (self, closest)
    for (i = lowerBound; i <= upperBound; i++) {
        double ulj = 0;
        if (ion[i].posvec.z < -0.5 * box.lz +
                              ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk. replacing 1 by ion[i].diameter -Yufei -Vikram -Vikram
        {
            PARTICLE dummy = PARTICLE(0, 0, 0, 0, 0, box.eout,
                                      VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, -0.5 * box.lz),
                                      box.lx, box.ly, box.lz);
            VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
            double r2 = r_vec.GetMagnitudeSquared();
            double d = 0.5 * ion[i].diameter;
            double d2 = d * d;
            double elj = 1.0;
            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double d6 = d2 * d2 * d2;
                ulj = 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
            }
        }
        lj_ion_leftdummy[i-lowerBound] = ulj;
    }

    // ion interacting with discretized left wall
//#pragma omp parallel default(shared) private(i, k, ulj, r2, r6, d, d2, d6, elj, wall_dummy, r_vec)
//    {
//#pragma omp for schedule(dynamic) nowait
//       for (i = lowerBound; i <= upperBound; i++) {
//            ulj = 0.0;
//            if (ion[i].posvec.z < -0.5 * box.lz +
//                                  ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk
//            {
//                for (k = 0; k < box.leftplane.size(); k++) {
//                    wall_dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.eout,
//                                          VECTOR3D(box.leftplane[k].posvec.x, box.leftplane[k].posvec.y,
//                                                   box.leftplane[k].posvec.z - 0.5 * ion[i].diameter), box.lx, box.ly,
//                                          box.lz);
//                    r_vec = ion[i].posvec - wall_dummy.posvec;
//
//                    if (r_vec.x > box.lx / 2) r_vec.x -= box.lx;
//                    if (r_vec.x < -box.lx / 2) r_vec.x += box.lx;
//                    if (r_vec.y > box.ly / 2) r_vec.y -= box.ly;
//                    if (r_vec.y < -box.ly / 2) r_vec.y += box.ly;
//
//                   r2 = r_vec.GetMagnitudeSquared();
//                    d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
//                    d2 = d * d;
//                    elj = 1.0;
//                    if (r2 < dcut2 * d2) {
//                        r6 = r2 * r2 * r2;
//                        d6 = d2 * d2 * d2;
//                        ulj += 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
//                    }
//                }
//            }
//            lj_ion_leftwall[i-lowerBound] = ulj;
//        }
//    }

    // right wall

    // ion interacting with right wall directly (self, closest)
    for (i = lowerBound; i <= upperBound; i++) {
        double ulj = 0;
        if (ion[i].posvec.z >
            0.5 * box.lz - ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk
        {
            PARTICLE dummy = PARTICLE(0, 0, 0, 0, 0, box.eout,
                                      VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, 0.5 * box.lz),
                                      box.lx, box.ly, box.lz);
            VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
            double r2 = r_vec.GetMagnitudeSquared();
            double d = 0.5 * ion[i].diameter;
            double d2 = d * d;
            double elj = 1.0;
            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double d6 = d2 * d2 * d2;
                ulj = 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
            }
        }
        lj_ion_rightdummy[i-lowerBound] = ulj;
    }

    // ion interacting with discretized right wall
//#pragma omp parallel default(shared) private(i, k, ulj, r2, r6, d, d2, d6, elj, wall_dummy, r_vec)
//    {
//#pragma omp for schedule(dynamic) nowait
//        for (i = lowerBound; i <= upperBound; i++) {
//            ulj = 0;
//            if (ion[i].posvec.z > 0.5 * box.lz -
//                                  ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk
//            {
//                for (k = 0; k < box.rightplane.size(); k++) {
//                    wall_dummy = PARTICLE(0, ion[i].diameter, 0, 0, 0, box.eout,
//                                          VECTOR3D(box.rightplane[k].posvec.x, box.rightplane[k].posvec.y,
//                                                   box.rightplane[k].posvec.z + 0.5 * ion[i].diameter), box.lx, box.ly,
//                                          box.lz);
//                    r_vec = ion[i].posvec - wall_dummy.posvec;
//
//                    if (r_vec.x > box.lx / 2) r_vec.x -= box.lx;
//                    if (r_vec.x < -box.lx / 2) r_vec.x += box.lx;
//                    if (r_vec.y > box.ly / 2) r_vec.y -= box.ly;
//                    if (r_vec.y < -box.ly / 2) r_vec.y += box.ly;
//
//                    r2 = r_vec.GetMagnitudeSquared();
//                    d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
//                    d2 = d * d;
//                    elj = 1.0;
//                    if (r2 < dcut2 * d2) {
//                        r6 = r2 * r2 * r2;
//                        d6 = d2 * d2 * d2;
//                        ulj += 4 * elj * (d6 / r6) * ((d6 / r6) - 1) + elj;
//                    }
//                }
//            }
//            lj_ion_rightwall[i-lowerBound] = ulj;
//        }
//    }
//coulomb interaction ion-rightwall///////////////////////////////////////////////////////////////////////////
#pragma omp parallel default(shared) private(i, k, wall_dummy, fqq_rightwall, fqq_csh_rightwall, dz_rightwall , r1_rightwall, r2_rightwall, fcsh_inf_rightwall, fcsh_z_rightwall, temp_vec_rightwall)
{
#pragma omp for schedule(dynamic) nowait
  for (i = lowerBound; i <= upperBound; i++)
  {
    fqq_rightwall = 0;
    fqq_csh_rightwall = 0;
    for (k = 0; k < box.rightplane.size(); k++) {
      wall_dummy = PARTICLE(0, 0, valency_counterion * -1, charge_meshpoint * 1.0, 0, box.eout, VECTOR3D(box.rightplane[k].posvec.x, box.rightplane[k].posvec.y, box.rightplane[k].posvec.z), box.lx, box.ly, box.lz);
      dz_rightwall  = ion[i].posvec.z - wall_dummy.posvec.z;
      r1_rightwall = sqrt(0.5 + (dz_rightwall  / box.lx) * (dz_rightwall  / box.lx));
      r2_rightwall = sqrt(0.25 + (dz_rightwall  / box.lx) * (dz_rightwall  / box.lx));
      fcsh_z_rightwall = 4 * box.lx * log((0.5 + r1_rightwall) / r2_rightwall) - fabs(dz_rightwall) * (2 * pi - 4 * atan(4 * fabs(dz_rightwall ) * r1_rightwall / box.lx));
      fcsh_inf_rightwall = -2 * pi * fabs(dz_rightwall);
      fqq_csh_rightwall += ion[i].q * (wall_dummy.q / (box.lx * box.lx)) * 0.5 * (1 / ion[i].epsilon + 1 / wall_dummy.epsilon) * (fcsh_inf_rightwall - fcsh_z_rightwall);

      temp_vec_rightwall = ion[i].posvec - wall_dummy.posvec;
      if (temp_vec_rightwall.x > box.lx / 2) temp_vec_rightwall.x -= box.lx;
      if (temp_vec_rightwall.x < -box.lx / 2) temp_vec_rightwall.x += box.lx;
      if (temp_vec_rightwall.y > box.ly / 2) temp_vec_rightwall.y -= box.ly;
      if (temp_vec_rightwall.y < -box.ly / 2) temp_vec_rightwall.y += box.ly;

      fqq_rightwall += ion[i].q * wall_dummy.q * 0.5 * (1.0 / ion[i].epsilon + 1.0 / wall_dummy.epsilon) / ((temp_vec_rightwall).GetMagnitude());

    }
    coulomb_rightwall[i-lowerBound] = fqq_rightwall + fqq_csh_rightwall;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//coulomb interaction ion-leftwall///////////////////////////////////////////////////////////////////////////
#pragma omp parallel default(shared) private(i, m, wall_dummy, fqq_leftwall, fqq_csh_leftwall, dz_leftwall, r1_leftwall, r2_leftwall, fcsh_inf_leftwall, fcsh_z_leftwall, temp_vec_leftwall)
{
#pragma omp for schedule(dynamic) nowait
    for (i = lowerBound; i <= upperBound; i++)
    {
      fqq_leftwall = 0;
      fqq_csh_leftwall = 0;
      for (m = 0; m < box.leftplane.size(); m++)
      {
        wall_dummy = PARTICLE(0, 0, valency_counterion * -1, charge_meshpoint * 1.0, 0, box.eout, VECTOR3D(box.leftplane[m].posvec.x, box.leftplane[m].posvec.y, box.leftplane[m].posvec.z), box.lx, box.ly, box.lz);

        dz_leftwall = ion[i].posvec.z - wall_dummy.posvec.z;
        r1_leftwall = sqrt(0.5 + (dz_leftwall / box.lx) * (dz_leftwall / box.lx));
        r2_leftwall = sqrt(0.25 + (dz_leftwall / box.lx) * (dz_leftwall / box.lx));
        fcsh_z_leftwall = 4 * box.lx * log((0.5 + r1_leftwall) / r2_leftwall) - fabs(dz_leftwall) * (2 * pi - 4 * atan(4 * fabs(dz_leftwall) * r1_leftwall / box.lx));
        fcsh_inf_leftwall = -2 * pi * fabs(dz_leftwall);
        fqq_csh_leftwall += ion[i].q * (wall_dummy.q / (box.lx * box.lx)) * 0.5 * (1 / ion[i].epsilon + 1 / wall_dummy.epsilon) * (fcsh_inf_leftwall - fcsh_z_leftwall);

        temp_vec_leftwall = ion[i].posvec - wall_dummy.posvec;
        if (temp_vec_leftwall.x > box.lx / 2) temp_vec_leftwall.x -= box.lx;
        if (temp_vec_leftwall.x < -box.lx / 2) temp_vec_leftwall.x += box.lx;
        if (temp_vec_leftwall.y > box.ly / 2) temp_vec_leftwall.y -= box.ly;
        if (temp_vec_leftwall.y < -box.ly / 2) temp_vec_leftwall.y += box.ly;
        fqq_leftwall += ion[i].q * wall_dummy.q * 0.5 * (1.0 / ion[i].epsilon + 1.0 / wall_dummy.epsilon) / ((temp_vec_leftwall).GetMagnitude());
      }
      coulomb_leftwall[i-lowerBound] = fqq_leftwall + fqq_csh_leftwall;
    }
  }

  // electrostatic interacting between right and left walls

  double total_coulomb_between_walls = 0;
  for (l = 0; l < box.rightplane.size(); l++)
  {
    double fqq_walls = 0.0;
    PARTICLE right_wall = PARTICLE(0, 0, valency_counterion * -1, charge_meshpoint * 1.0, 0, box.eout, VECTOR3D(box.rightplane[l].posvec.x, box.rightplane[l].posvec.y, box.rightplane[l].posvec.z), box.lx, box.ly, box.lz);
    for (b = 0; b < box.leftplane.size(); b++)
    {
        if (l == b) continue;
      PARTICLE left_wall = PARTICLE(0, 0, valency_counterion * -1, charge_meshpoint * 1.0, 0, box.eout, VECTOR3D(box.leftplane[b].posvec.x, box.leftplane[b].posvec.y, box.leftplane[b].posvec.z), box.lx, box.ly, box.lz);

      VECTOR3D r_vec_walls = right_wall.posvec - left_wall.posvec;
      fqq_walls = fqq_walls + 0.5 * right_wall.q * left_wall.q * 0.5 * (1.0 / right_wall.epsilon + 1.0 / left_wall.epsilon) / ((r_vec_walls).GetMagnitude());
  }
  total_coulomb_between_walls += fqq_walls;
}

    double ion_ion = 0;
    double total_lj_ion_ion = 0;
    double total_lj_ion_leftdummy = 0;
  //  double total_lj_ion_leftwall = 0;
    double total_lj_ion_rightdummy = 0;
  //  double total_lj_ion_rightwall = 0;
    double total_coulomb_leftwall = 0;
    double total_coulomb_rightwall = 0;
    for (i = lowerBound; i <= upperBound; i++) {
        ion_ion += ion_energy[i-lowerBound];
        total_lj_ion_ion += lj_ion_ion[i-lowerBound];
        total_lj_ion_leftdummy += lj_ion_leftdummy[i-lowerBound];
//        total_lj_ion_leftwall += lj_ion_leftwall[i-lowerBound];
        total_lj_ion_rightdummy += lj_ion_rightdummy[i-lowerBound];
//        total_lj_ion_rightwall += lj_ion_rightwall[i-lowerBound];
        total_coulomb_rightwall += coulomb_rightwall[i-lowerBound];
        total_coulomb_leftwall += coulomb_leftwall[i-lowerBound];
    }
    // electrostatic potential energy
    double coulomb = (ion_ion) * scalefactor;    // scalefactor is there to ensure proper units; defined in utility.h
    total_coulomb_rightwall = (total_coulomb_rightwall) * scalefactor;
    total_coulomb_leftwall = (total_coulomb_leftwall) * scalefactor;
    total_coulomb_between_walls = (total_coulomb_between_walls) * scalefactor;

    total_lj_ion_ion = 0.5 * total_lj_ion_ion;
    // factor of half for double counting, same reasoning as electrostatic energy

    potential = coulomb + total_lj_ion_ion + total_lj_ion_leftdummy + total_lj_ion_rightdummy + total_coulomb_rightwall + total_coulomb_leftwall;
  //  potential = coulomb + total_lj_ion_ion + total_lj_ion_leftwall + total_lj_ion_rightwall;

    //MPI Operations
    if (world.size() > 1) {
        //broadcasting using all_reduce = reduce + broadcast
        //void all_reduce(const communicator & comm, const T & in_value, T & out_value, Op op);
        all_reduce(world, potential, totalPotential,std::plus<double>());
    }else{

        totalPotential=potential;
    }

    totalPotential = totalPotential + total_coulomb_between_walls;

    return totalPotential;
}
