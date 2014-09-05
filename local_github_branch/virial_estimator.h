//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 4 Sep, 2014
//
//Virial Estimator is a class which calculates the average and instantenous
//energies of the system using the differences of the potential felt by each
//of the "cyclic paths" (e.g., the locations of the beads on the necklace
//with respect to the central potential). 
//
//This is an energy estimator SPECIFICALLY DESIGNED FOR THE HARMONIC OSCILLATOR

#ifndef VIRIAL_ESTIMATOR_H
#define VIRIAL_ESTIMATOR_H

#include "necklace.h"

class virial_estimator{

public:
    void initialize_estimator(int number_of_beads, long double spring_constant);
    void update_energy(necklace my_necklace, int weight_of_average);
    long double get_average_energy(){return avg_estimated_energy_;}
    long double get_instant_energy(){return instantaneous_estimated_energy_;}

private:
    long double avg_estimated_energy_;
    long double instantaneous_estimated_energy_;
    long double multiplier_;
    long double spring_constant_;
};

#endif
