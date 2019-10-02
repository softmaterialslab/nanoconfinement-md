// This is vertex class

#ifndef _VERTEX_H
#define _VERTEX_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"

// NOTE: Some members and member functions of the VERTEX class are not relevant to the system that this code is designed to simulate.
// But they are retained to help readability of the CPMD code where they become critical.

class VERTEX
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & posvec;
      ar & a;
      ar & normalvec;
      ar & r;
      ar & theta;
      ar & phi;
      ar & mu;
      ar & w;
      ar & vw;
      ar & fw;
      ar & ke;
      ar & pe;
      ar & wmean;
      ar & presumgwEw;
      ar & presumgEwEq;
      ar & presumgEwEw;
      ar & presumfwEw;
      ar & presumfEwEq;
      ar & presumhEqEw;
      ar & Greens;
      ar & ndotGradGreens;
      ar & Gion;
      ar & gradGion;
      ar & q;
      ar & epsilon;
    }

  public:

  // members
  VECTOR3D posvec;			// position vector of the vertex
  long double a;			// area of the vertex
  VECTOR3D normalvec;     		// normal vector on the surface pointing interior to exterior
  double r, theta, phi;			// polar coordinates of the vertex
  long double mu;			// 'mass' of the induced charge
  long double w;			// discretized induced charge on the vertex
  long double vw;			// 'velocity' of the induced charge
  long double fw; 			// 'force' on the induced charge
  long double ke;              		// 'kinetic energy' of the induced charge
  double pe;				// potential energy of the fake degrees associated with vertices
  double wmean;				// mean w computed on fmd
  double q;		// charge of the meshpoint
  double epsilon;	// dielectric constant of the medium

  // member vectors
  vector<long double> presumgwEw;	// result of a precalculation-- gwEw
  vector<long double> presumgEwEq;	// result of a precalculation-- gEwEq
  vector<long double> presumgEwEw;	// result of a precalculation-- gEwEw
  vector<long double> presumfwEw;	// result of a precalculation-- fwEw
  vector<long double> presumfEwEq;	// result of a precalculation-- fwEw
  vector<long double> presumhEqEw;	// result of a precalculation-- hEqEw
  vector<long double> Greens;		// precalculate Gww
  vector<long double> ndotGradGreens;	// precalculate n.gradGww
  vector<long double> Gion;		// Greens function between induced charge and ion
  vector<VECTOR3D> gradGion;		// gradient of the above

  // member functions

  // make a vertex
  VERTEX(VECTOR3D initial_position = VECTOR3D(0,0,0), double initial_charge = 0, double initial_diconst = 0, double initial_area = 0, VECTOR3D initial_normal = VECTOR3D(0,0,0))
  {
    posvec = initial_position;
    q = initial_charge;
    epsilon = initial_diconst;
    a = initial_area;
    normalvec = initial_normal;
  }

  // update position of fake degree of freedom (auxillary variable, fake variable)
  void update_position(double dt)
  {
    w = w + dt * vw;
    return;
  }

  // update velocity of fake degree of freedom (in simple velocity verlet; no thermostat)
  void update_velocity(double dt)
  {
    vw = vw + 0.5 * dt * fw / (mu);
    return;
  }

  // update velocity with integrator that unifies velocity verlet and Nose-Hoover
  void new_update_velocity(double dt, THERMOSTAT main_bath, long double expfac)
  {
    vw = vw * expfac + 0.5 * dt * (fw / (mu)) * sqrt(expfac);
    return;
  }

  // compute kinetic energy of fake degree
  void kinetic_energy()
  {
    ke = 0.5 * mu * vw * vw;
    return;
  }

  // get the polar coordinates for the vertex
  void get_polar()
  {
    r = posvec.GetMagnitude();
    theta = posvec.y > 0 ? acos(posvec.z / r) : 2 * pi - acos(posvec.z / r);
    phi = posvec.y > 0 ? acos(posvec.x / r) : - acos(posvec.x / r);
    return;
  }
};

#endif
