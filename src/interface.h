// This is a header file for the INTERFACE class.

#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "utility.h"
#include "vector3d.h"
#include "particle.h"
#include "vertex.h"

class INTERFACE
{
  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & posvec;
      ar & ein;
      ar & eout;
      ar & em;
      ar & ed;
      ar & lx;
      ar & ly;
      ar & lz;
      ar & lB_in;
      ar & lB_out;
      ar & inv_kappa_in;
      ar & inv_kappa_out;
      ar & mean_sep_in;
      ar & mean_sep_out;
      ar & width;
      ar & leftplane;
      ar & rightplane;

    }

  public:

  VECTOR3D posvec;		// position vector of the inteface, for sphere its center
  double ein; 			// permittivity of inside medium
  double eout; 			// permittivity of outside medium
  double em;			// permittivity at the interface, mean
  double ed;			// permittivity change at the inteface, difference scaled by 4 pi
  double lx;			// length of the box in x direction
  double ly;			// length of the box in y direction
  double lz; 			// length of the box in z direction
  double lB_in;			// Bjerrum length inside
  double lB_out;		// Bjerrum length outside
  double inv_kappa_in;		// debye length inside
  double inv_kappa_out;		// debye length outside
  double mean_sep_in;		// mean separation inside
  double mean_sep_out;		// mean separation outside

  double width;			// discretization width
  vector<VERTEX> leftplane;
  vector<VERTEX> rightplane;

  void set_up(double, double, int, int, double, double, double);
  void put_saltions_inside(vector<PARTICLE>&, int, int, double, double, double, vector<PARTICLE>&, unsigned int, int, double, double);
  void discretize(double, double);
  void generate_lammps_datafile(vector<PARTICLE>&, int, int, vector<PARTICLE>&, double, double, int, int, double, double);

  INTERFACE(VECTOR3D initial_posvec = VECTOR3D(0,0,0), double initial_ein = 1, double initial_eout = 1)
  {
    posvec = initial_posvec;
    ein = initial_ein;
    eout = initial_eout;
  }

  // total charge inside
  double total_charge_inside(vector<PARTICLE>& ion)
  {
    double charge = 0;
    for (unsigned long i = 0; i < ion.size(); i++)
    {
      double r = 0.5 * ion[i].diameter;
      if ((ion[i].posvec.x <= lx - r && ion[i].posvec.x >= -lx + r) &&
          (ion[i].posvec.y <= ly - r && ion[i].posvec.y >= -ly + r) &&
          (ion[i].posvec.z <= lz - r && ion[i].posvec.z >= -lz + r))
        charge += ion[i].q;
      else
      {
        charge = -1;
        break;
      }
    }
    return charge;
  }

  // total induced charge (should be 0 for this system)
  // this is not relevant in uniform dielectric but is retained to refer back to when looking at cpmd code for non-uniform media
  double total_induced_charge(vector<VERTEX>& s)
  {
    double charge = 0;
    for (unsigned int k = 0; k < s.size(); k++)
      charge += s[k].w * s[k].a;
    return charge;
  }

};

#endif
