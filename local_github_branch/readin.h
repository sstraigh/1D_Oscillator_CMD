//Author: Shelby Straight, Paesani lab
//Date Last Edited: 4, Sep 2014
//
//readin.h/.cpp define the object responsible for reading and setting 
//experimental parameters (hbar, physical mass, spring constant, etc..)
//from the file "parameters.dat"

#ifndef READIN_H
#define READIN_H

class readin{

public:
void read_parameters(int& Number_of_Steps, int& Number_of_Beads,
		     int& Thermostat_Length, long double& Temperature,
		     int& Capture_Rate, long double& Gamma, int& Potential_Choice,
		     int& Detach_Centroid_Thermostat, int& Continue_Previous_Simulation,
		     long double& Hbar, long double& Central_Spring_Constant,
		     long double& Particle_Mass, long double& Initial_Velocity,
		     long double& Initial_Position, long double& Timestep,
		     int& Spatial_Deg_of_Freedom); 

};

#endif
