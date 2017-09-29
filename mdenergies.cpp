// This file contains the routine to compute the total potential energy of the system
// Electrostatic and Excluded volume (LJ) contributions

#include "energies.h"

// Potential energy
long double energy_functional(vector<PARTICLE>& ion, INTERFACE& box)
{
  unsigned int i, j, k;
  double fqq, fqq_csh;
  double fcsh_z, r1, r2, dz, fcsh_inf;
  VECTOR3D temp_vec;
 
  vector<double> ion_energy(ion.size(), 0.0);
  // push_back is not compatible with pragma, so we have initialize vector in this way.
  
  // charged sheets method is used to compute Coulomb interactions; an Ewald version should be designed to compare and ensure that long-range effects are taken into account in either methods
  
#pragma omp parallel default(shared) private(i, j, fqq, fqq_csh, dz, r1, r2, fcsh_inf, fcsh_z,temp_vec)
{
  #pragma omp for schedule(dynamic) nowait
  for (i = 0; i < ion.size(); i++)
  {
    fqq = 0;
    fqq_csh =0;
    for (j = 0; j < ion.size(); j++)
    {
      dz = ion[i].posvec.z-ion[j].posvec.z;
      r1 = sqrt(0.5 + (dz/box.lx)*(dz/box.lx));
      r2 = sqrt(0.25 + (dz/box.lx)*(dz/box.lx));
      fcsh_z = 4 * box.lx * log((0.5 + r1)/r2) - fabs(dz) * (2 * pi - 4 * atan(4*fabs(dz)*r1/box.lx));
      fcsh_inf = -2 * pi * fabs(dz);
      fqq_csh += ion[i].q * (ion[j].q/ (box.lx * box.lx)) * 0.5 * (1/ion[i].epsilon + 1/ion[j].epsilon) * (fcsh_inf - fcsh_z);

      if (i == j) continue;

      temp_vec=ion[i].posvec-ion[j].posvec;

      if (temp_vec.x > box.lx/2) temp_vec.x -= box.lx;
      if (temp_vec.x < -box.lx/2) temp_vec.x += box.lx;
      if (temp_vec.y > box.ly/2) temp_vec.y -= box.ly;
      if (temp_vec.y < -box.ly/2) temp_vec.y += box.ly;

      fqq += 0.5 * ion[i].q * ion[j].q * 0.5 * (1.0 / ion[i].epsilon + 1.0 / ion[j].epsilon) / ((temp_vec).GetMagnitude());
    }
    ion_energy[i] = fqq + fqq_csh;
  }
}

  double ion_ion = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    ion_ion = ion_ion + ion_energy[i];
  
  // electrostatic potential energy
  double coulomb = (ion_ion) * scalefactor;	// scalefactor is there to ensure proper units; defined in utility.h

  // Excluded volume interaction energy given by purely repulsive LJ
  
  // ion-ion
  vector<double> lj_ion_ion(ion.size(), 0.0);
  VECTOR3D r_vec;
  PARTICLE wall_dummy;
  double ulj, uljcc, elj, r, r6, d, d2, d6;
  #pragma omp parallel default(shared) private(i, j, uljcc, r, r2, r6, d, d2, d6, elj, r_vec)
  {
    #pragma omp for schedule(dynamic) nowait
  for (i = 0; i < ion.size(); i++)
  {
    uljcc = 0.0;
    for (j = 0; j < ion.size(); j++)
    {
      if (j == i) continue;
      r_vec = ion[i].posvec - ion[j].posvec;
      
      if (r_vec.x > box.lx/2) r_vec.x -= box.lx;
      if (r_vec.x < -box.lx/2) r_vec.x += box.lx;
      if (r_vec.y > box.ly/2) r_vec.y -= box.ly;
      if (r_vec.y < -box.ly/2) r_vec.y += box.ly; 
      
      r = r_vec.GetMagnitude();
      d = 0.5 * (ion[i].diameter + ion[j].diameter);
      elj = 1.0;
      if (r < dcut * d)
      {
	r2 = r * r;
	r6 = r2 * r2 * r2;
	d2 = d * d;
	d6 = d2 * d2 * d2;
	uljcc = uljcc +  4 * elj * (d6 / r6) * ( ( d6 / r6 ) - 1 ) + elj;
      }
      else
	uljcc = uljcc + 0.0;
    }
    lj_ion_ion[i] =uljcc;
  }
  }
  // ion-box
  
  // left wall
  
  // ion interacting with left wall directly (self, closest)
  vector<double> lj_ion_leftdummy(ion.size(), 0.0);
  for (unsigned int i = 0; i < ion.size(); i++)
  {
    double ulj = 0;
    if(ion[i].posvec.z < -0.5*box.lz + ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk. replacing 1 by ion[i].diameter -Yufei -Vikram -Vikram
    {
      PARTICLE dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.eout,VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, -0.5*box.lz - 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
      VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
      double r2 = r_vec.GetMagnitudeSquared();
      double d = 0.5 * (ion[i].diameter + dummy.diameter);
      double d2 = d * d;
      double elj = 1.0;
      if (r2 < dcut2 * d2)
      {
	double r6 = r2 * r2 * r2;
	double d6 = d2 * d2 * d2;
	ulj = 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
      }
    }
    lj_ion_leftdummy[i] = ulj;
  }
  
  // ion interacting with discretized left wall
  vector<double> lj_ion_leftwall(ion.size(), 0.0);
  #pragma omp parallel default(shared) private(i, k, ulj, r2, r6, d, d2, d6, elj, wall_dummy,r_vec)
  {
    #pragma omp for schedule(dynamic) nowait
  for (i = 0; i < ion.size(); i++)
  {
    ulj = 0.0;
    if(ion[i].posvec.z < -0.5*box.lz + ion[i].diameter)   // avoiding calculating interactions between left wall and ions in bulk
    {
      for (k = 0; k < box.leftplane.size(); k++)
      {
	wall_dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.eout, VECTOR3D(box.leftplane[k].posvec.x, box.leftplane[k].posvec.y, box.leftplane[k].posvec.z - 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
	r_vec = ion[i].posvec - wall_dummy.posvec;
	
	if (r_vec.x > box.lx/2) r_vec.x -= box.lx;
	if (r_vec.x < -box.lx/2) r_vec.x += box.lx;
	if (r_vec.y > box.ly/2) r_vec.y -= box.ly;
	if (r_vec.y < -box.ly/2) r_vec.y += box.ly;
	
	r2 = r_vec.GetMagnitudeSquared();
	d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
	d2 = d * d; 
	elj = 1.0;
	if (r2 < dcut2 * d2)
	{
	  r6 = r2 * r2 * r2;
	  d6 = d2 * d2 * d2;
	  ulj += 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
	}
      }
    }
    lj_ion_leftwall[i] = ulj;
  }
  }
  
  // right wall
  
  // ion interacting with right wall directly (self, closest)
  vector<double> lj_ion_rightdummy(ion.size(), 0.0);
  for (unsigned int i = 0; i < ion.size(); i++)
  {
    double ulj = 0;
    if(ion[i].posvec.z > 0.5*box.lz - ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk
    {
      PARTICLE dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.eout,VECTOR3D(ion[i].posvec.x, ion[i].posvec.y, 0.5*box.lz + 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
      VECTOR3D r_vec = ion[i].posvec - dummy.posvec;
      double r2 = r_vec.GetMagnitudeSquared();
      double d = 0.5 * (ion[i].diameter + dummy.diameter);
      double d2 = d * d;
      double elj = 1.0;
      if (r2 < dcut2 * d2)
      {
	double r6 = r2 * r2 * r2;
	double d6 = d2 * d2 * d2;
	ulj = 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
      }
    }
    lj_ion_rightdummy[i] = ulj;
  }
  
  // ion interacting with discretized right wall
  vector<double> lj_ion_rightwall(ion.size(), 0.0);
  #pragma omp parallel default(shared) private(i, k, ulj, r2, r6, d, d2, d6, elj, wall_dummy,r_vec)
  {
    #pragma omp for schedule(dynamic) nowait
  for (i = 0; i < ion.size(); i++)
  {
    ulj = 0;
    if(ion[i].posvec.z > 0.5*box.lz - ion[i].diameter)  // avoiding calculating interactions between right wall and ions in bulk
    {
      for (k = 0; k < box.rightplane.size(); k++)
      {
	wall_dummy = PARTICLE(0,ion[i].diameter,0,0,0,box.eout, VECTOR3D(box.rightplane[k].posvec.x, box.rightplane[k].posvec.y, box.rightplane[k].posvec.z + 0.5*ion[i].diameter),box.lx,box.ly,box.lz);
	r_vec = ion[i].posvec - wall_dummy.posvec;
	
	if (r_vec.x > box.lx/2) r_vec.x -= box.lx;
	if (r_vec.x < -box.lx/2) r_vec.x += box.lx;
	if (r_vec.y > box.ly/2) r_vec.y -= box.ly;
	if (r_vec.y < -box.ly/2) r_vec.y += box.ly;
	
	r2 = r_vec.GetMagnitudeSquared();
	d = 0.5 * (ion[i].diameter + wall_dummy.diameter);
	d2 = d * d;
	elj = 1.0;
	if (r2 < dcut2 * d2)
	{
	  r6 = r2 * r2 * r2;
	  d6 = d2 * d2 * d2;
	  ulj += 4 * elj * (d6 / r6) * ( (d6 / r6) - 1 ) + elj;
	}
      }
    }
    lj_ion_rightwall[i] = ulj;
  }
  }
  
  double total_lj_ion_ion = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_ion += lj_ion_ion[i];
  total_lj_ion_ion = 0.5 * total_lj_ion_ion;
  // factor of half for double counting, same reasoning as electrostatic energy 
  
  double total_lj_ion_leftdummy = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_leftdummy += lj_ion_leftdummy[i];
  
  double total_lj_ion_leftwall = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_leftwall += lj_ion_leftwall[i];
  
  double total_lj_ion_rightdummy = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_rightdummy += lj_ion_rightdummy[i];
  
  double total_lj_ion_rightwall = 0;
  for (unsigned int i = 0; i < ion.size(); i++)
    total_lj_ion_rightwall += lj_ion_rightwall[i];

  double potential = coulomb + total_lj_ion_ion + total_lj_ion_leftdummy + total_lj_ion_leftwall 
			    + total_lj_ion_rightdummy + total_lj_ion_rightwall;
		      
  return potential;
}

