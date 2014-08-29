#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "thermostat.h"

class particle{
    
public:
    particle();
    
    long double kinetic_energy;
    long double potential_energy;
    long double location; //could be an array of three values, but because
		            //we are dealing with the special 1D case it is simply one var
			    //(3 value array would also require the velocity to be expressed
			    //as a vector rather than a scalar)

    long double Beta;    
    long double omega_k;
    long double mass;
    long double velocity;

    thermostat nhc;
    long double acceleration;

    void initialize(long double mas, long double displacement, long double spr_constant, long double velocity, int ndof, int chain_length, long double timestep, long double kT, long double spring_freq);

    void integrate_thermostat(long double timestep);

    void deallocate();
};

#endif
