//Author: Shelby Straight, Paesani Lab
//Date last edited: Sep 3, 2014
//
//Gaussian.h is the header file for a "gaussian" object,
//which is instantiated when a some part of the simulaiton
//(i.e., sampling of initial velocities with respect to
//kT/sqrt(mass)) requires a random number generator

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

class gaussian{

public:
    long double compute_gaussian(long double vel);

};
#endif
