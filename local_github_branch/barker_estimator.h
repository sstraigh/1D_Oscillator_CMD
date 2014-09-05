//Author: Shelby Straight, Paesani Lab
//Date Last Edited: 4 Sep, 2014
//
//Barker Estimator is a class which calculates the average and instantenous
//energies of the system using the differences of "cyclic paths" (e.e., the
//difference in location of the beads on the necklace with respect to their
//nearest neighbors). 
//
//This is an energy estimator SPECIFICALLY DESIGNED FOR THE HARMONIC OSCILLATOR
//
//It is different than the virial estimator in the following ways:
//     -it uses the differences in location of the beads rather than the
//         potential felt by individual beads to estimate the energy
//     -it is not as accurate as the virial estimator, and its error increases
//         with P, the number of beads (Tuckerman, Stat Mech)

#ifndef BARKER_ESTIMATOR_H
#define BARKER_ESTIMATOR_H

#include "necklace.h"

class barker_estimator{

public:
    long double get_avg_KE(){return avg_kinetic_energy_;}
    long double get_avg_PE(){return avg_potential_energy_;}
    long double get_instant_PE(){return instant_potential_energy_;}
    long double get_instant_KE(){return instant_kinetic_energy_;}
    long double get_instant_total_E(){return (instant_potential_energy_+instant_kinetic_energy_);}
    long double get_avg_total_E(){return (avg_kinetic_energy_+avg_potential_energy_);}

    void initialize_estimator(int number_of_beads, long double spring_constant,
			      long double particle_mass, long double kT, 
			      long double hbar);

    void update_energy(necklace my_necklace, int weight_of_average);

private:
    long double avg_kinetic_energy_;
    long double avg_potential_energy_;
    long double instant_kinetic_energy_;
    long double instant_potential_energy_;
    long double kinetic_prefactor_;
    long double spring_constant_;

};

#endif
