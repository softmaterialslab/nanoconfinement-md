#ifndef _FORCES_H
#define _FORCES_H

#include "vertex.h"
#include "particle.h"
#include "interface.h"
#include "functions.h"

void for_md_calculate_force(vector<PARTICLE>&, INTERFACE&, char);

#endif 
