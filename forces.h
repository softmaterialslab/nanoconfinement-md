// This is a header file that contains functions that compute the forces on particles

#ifndef _FORCES_H
#define _FORCES_H

#include "vertex.h"
#include "particle.h"
#include "interface.h"
#include "functions.h"

void for_md_calculate_force(vector<PARTICLE>&, INTERFACE&, char, vector<vector<VECTOR3D> >&, vector <VECTOR3D> &, vector <VECTOR3D> &, vector <VECTOR3D> &, vector <VECTOR3D> &, vector <VECTOR3D> &, vector <VECTOR3D> &);

#endif 
