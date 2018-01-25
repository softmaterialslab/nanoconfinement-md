// This is particle class

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"

class PARTICLE 
{
  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & id;
      ar & diameter;
      ar & valency;
      ar & q;
      ar & m;
      ar & epsilon;
      ar & posvec;
      ar & velvec;
      ar & forvec;
      ar & pe;
      ar & ke;
      ar & energy;
      ar & lx;
      ar & ly;
      ar & lz;

    }

  public:

  // members
  int id;		// id of the particle
  double diameter;	// diameter of the particle
  int valency;		// valency of the ion
  double q;		// charge of the particle
  double m; 		// mass of the particle
  double epsilon;	// dielectric constant of the medium
  VECTOR3D posvec;	// position vector of the particle
  VECTOR3D velvec;	// velocity vector of the particle
  VECTOR3D forvec;	// force vector on the particle
  double pe;		// potential energy
  long double ke;	// kinetic energy
  double energy;	// energy
  double lx;		// box length in x
  double ly;		// box length in y
  double lz;		// box length in z
  
  // member functions
  
  // make a particle
  PARTICLE(int initial_id = 0, double initial_diameter = 0, int initial_valency = 0, double initial_charge = 0, double initial_mass = 0, double initial_diconst = 0, VECTOR3D initial_position = VECTOR3D(0,0,0), double initial_lx = 0,double initial_ly = 0, double initial_lz = 0)
  {
    id = initial_id;
    diameter = initial_diameter;
    valency = initial_valency;
    q = initial_charge;
    m = initial_mass;
    epsilon = initial_diconst;
    posvec = initial_position;
    lx = initial_lx;
    ly = initial_ly;
    lz = initial_lz;
  }
  // other quantities in the particle class are not initialized when you make a particle. they should be initialized separately when utilized in functions.
  
  // update position of the particle
  void update_position(double dt)			
  {
    posvec = ( posvec + (velvec ^ dt) );
    // the following lines implement periodic boundary conditions in x and y direction; system confined in z
    // while is a better condition than if because it ensures that the particle comes back to the main cell; but it is costly
    if(posvec.x > lx/2)
      posvec.x -= lx;
    if(posvec.x < -lx/2)
      posvec.x += lx;
    if(posvec.y > ly/2)
      posvec.y -= ly;
    if(posvec.y < -ly/2)
      posvec.y += ly;
    return;
  }
  
  // update velocity of the particle (in simple velocity verlet; no thermostat)
  void update_velocity(double dt)	
  {
    velvec = ( velvec + ( forvec ^ ( 0.5 * dt ) ) );
    return;
  }
  
  // update velocity with integrator that unifies velocity verlet and Nose-Hoover
  void new_update_velocity(double dt, THERMOSTAT main_bath, long double expfac)
  {
    velvec = ( ( velvec ^ (expfac)  ) + ( forvec ^ (0.5 * dt * sqrt(expfac)) ) );
    return;
  }
  
  // calculate kinetic energy of a particle
  void kinetic_energy()				
  {
    ke = 0.5 * m * velvec.GetMagnitude() * velvec.GetMagnitude();
    return;
  }
};

#endif
