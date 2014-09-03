//Author: Shelby Straight, Paesani Lab
//Date Last Edited: Sep 3, 2014
//
//Particle.h is the header file which defines class particle-
//each particle can be thought of as a "bead" in a "necklace",
//and the necklace captures the quantum dynamics of the particle
//it is used to model.
//
//TODO: rename particle class to bead?


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
  long double get_location(){return location;}
  long double set_location(long double loc){location=loc;}

  long double get_velocity(){return velocity;}
  long double set_velocity(long double vel){velocity=vel;}

  long double get_acceleration(){return acceleration;}
  long double set_acceleration(long double acc){acceleration=acc;}

  long double get_mass(){return mass;}

  thermostat nhc;

  void initialize(long double mas, long double displacement, 
		  long double spr_constant, long double velocity, 
		  int ndof, int chain_length, long double timestep, 
		  long double kT, long double spring_freq);

  void integrate_thermostat(long double timestep);

  void deallocate();

private:
  long double velocity;
  long double location;
  long double acceleration;

  long double Beta;
  long double omega_k;
  long double mass;
};

#endif
