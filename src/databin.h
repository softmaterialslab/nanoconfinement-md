// This is bin class

#ifndef _DATABIN_H
#define _DATABIN_H

#include "utility.h"

class DATABIN
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & n;
        ar & width;
        ar & volume;
        ar & lower;
        ar & higher;
    }


public:
    
    // members
    int n;		// number of ions in the bin
    double width;	// width of the bin
    double volume;	// volume of the bin
    double lower;	// lower value of bin
    double higher;	// higher value of bin
    
    // member functions
    
    // make a bin
    DATABIN(int initial_n = 0, double initial_width = 0, double initial_volume = 0, double initial_lower = 0, double initial_higher = 0 )
    {
      n = initial_n;
      width = initial_width;
      volume = initial_volume;
      lower = initial_lower;
      higher = initial_higher;
    }
    
    // set up bin
    void set_up(int bin_num, double bin_width, double lx, double ly, double lz)
    {
      n = 0;			// initialize number of ions in bin to be zero
      width = bin_width;
      volume = (bin_width * lx * ly) * (unitlength * unitlength * unitlength) * 0.6022;
      lower = 0.5*(-lz) + bin_num * bin_width;
      higher = 0.5*(-lz) + (bin_num + 1) * bin_width;
    }
};

#endif
