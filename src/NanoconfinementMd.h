// This is a header file for the INTERFACE class.

#ifndef _NanoconfinementMd_H
#define _NanoconfinementMd_H

#include <boost/program_options.hpp>
#include "utility.h"
#include "interface.h"
#include "particle.h"
#include "vertex.h"
#include "databin.h"
#include "control.h"
#include "functions.h"
#include "thermostat.h"

void md(vector<PARTICLE>&, INTERFACE&, vector<THERMOSTAT>&, vector<DATABIN>&, CONTROL&, string&, double&);

using namespace boost::program_options;

class NanoconfinementMd
{
  public:
    int startSimulation(int argc, char *argv[], bool XY_Map);
};

#endif
