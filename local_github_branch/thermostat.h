//Author: Shelby Straight, Paesani Lab
//Date last edited: 3 Sep 2014
//
//Thermostat refers to the Nose Hoover Chain implementation of a "bath" which
//is coupled to each bead in the necklace-i.e., to each degree of freedom in
//the simulation.
//
//This has the effect of placing the particle to which it is coupled in the
//NVT (canonical) ensemble. Removing the thermostat replaces the particle in
//the NVE (microcanonical) ensemble. This is done in the last CMD simulation,
//where the thermostat is removed from the centroid.

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>

class thermostat{

 public:

  void initialize_thermostat(int ndof, int chain_length, long double timestep, 
			     long double kT, long double particle_velocity, 
			     long double particle_mass, 
			     long double spring_freq);

  void deallocate();

  long double calculate_Nhamiltonian();

  //accessed in classes which "have a" thermostat, therefore
  //set to public
  int chain_length; //M in Tuckerman(2010), "size" of the chain

  //the following arrays are public as they are used in the particle
  //class's integrate_thermostat() routine, which scales the particle
  //velocity according to the dynamics of the nose-hoover-chain
  long double * nose_accel;
  long double * nose_pos;
  long double * nose_velocity;
  long double * Q_param;
  long double suzuki_coefficients[5];

 private:
  long double p_velocity;
  long double p_mass;
  long double tau;
  long double Beta;

  int degrees_of_freedom;

  //contribution of the thermostat chain to the extended hamiltonian for
  //the system under consideration. Returned by public member function
  //calculate_Nhamiltonian()
  long double Nhamiltonian;

  //private variables corresponding to the positional and kinetic portions
  //of the hamiltonian for this particle
  long double N_pos_H;
  long double N_kin_H;

};

#endif
